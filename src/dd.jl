const ROSENBROCK_INV_CUTOFF = 7 # https://github.com/SciML/OrdinaryDiffEq.jl/pull/1539

struct StaticWOperator{isinv, T}
    W::T
    function StaticWOperator(W::T, callinv = true) where {T}
        isinv = size(W, 1) <= ROSENBROCK_INV_CUTOFF

        # when constructing W for the first time for the type
        # inv(W) can be singular
        _W = if isinv && callinv
            inv(W)
        else
            W
        end
        new{isinv, T}(_W)
    end
end
isinv(W::StaticWOperator{S}) where {S} = S
Base.:\(W::StaticWOperator, v) = isinv(W) ? W.W * v : W.W \ v
Base.:\(W::StaticWOperator, v::AbstractArray) = isinv(W) ? W.W * v : W.W \ v

function calc_tderivative!(integrator, cache, dtd1, repeat_step)
    @inbounds begin
        @unpack t, dt, uprev, u, f, p = integrator
        @unpack du2, fsalfirst, dT, tf, linsolve_tmp = cache

        # Time derivative
        if !repeat_step # skip calculation if step is repeated
            if DiffEqBase.has_tgrad(f)
                f.tgrad(dT, uprev, p, t)
            else
                tf.uprev = uprev
                tf.p = p
                derivative!(dT, tf, t, du2, integrator, cache.grad_config)
            end
        end

        @.. broadcast=false linsolve_tmp=fsalfirst + dtd1 * dT
    end
end

function calc_tderivative!(integrator::ODEIntegrator{algType, IIP, <:Array}, cache, dtd1,
                           repeat_step) where {algType, IIP}
    @inbounds begin
        @unpack t, dt, uprev, u, f, p = integrator
        @unpack du2, fsalfirst, dT, tf, linsolve_tmp = cache

        # Time derivative
        if !repeat_step # skip calculation if step is repeated
            if DiffEqBase.has_tgrad(f)
                f.tgrad(dT, uprev, p, t)
            else
                tf.uprev = uprev
                if !(p isa DiffEqBase.NullParameters)
                    tf.p = p
                end
                derivative!(dT, tf, t, du2, integrator, cache.grad_config)
            end
        end

        @inbounds @simd ivdep for i in eachindex(uprev)
            linsolve_tmp[i] = fsalfirst[i] + dtd1 * dT[i]
        end
    end
end

function calc_tderivative(integrator, cache)
    @unpack t, dt, uprev, u, f, p = integrator

    # Time derivative
    if DiffEqBase.has_tgrad(f)
        dT = f.tgrad(uprev, p, t)
    else
        tf = cache.tf
        tf.u = uprev
        tf.p = p
        dT = derivative(tf, t, integrator)
    end
    dT
end

"""
    calc_J(integrator, cache, next_step::Bool = false)

Return a new Jacobian object.

If `integrator.f` has a custom Jacobian update function, then it will be called. Otherwise,
either automatic or finite differencing will be used depending on the `uf` object of the
cache. If `next_step`, then it will evaluate the Jacobian at the next step.
"""
function calc_J(integrator, cache, next_step::Bool = false)
    @unpack dt, t, uprev, f, p, alg = integrator
    if next_step
        t = t + dt
        uprev = integrator.u
    end

    if alg isa DAEAlgorithm
        if DiffEqBase.has_jac(f)
            J = f.jac(duprev, uprev, p, t)
        else
            @unpack uf = cache
            x = zero(uprev)
            J = jacobian(uf, x, integrator)
        end
    else
        if DiffEqBase.has_jac(f)
            J = f.jac(uprev, p, t)
        else
            @unpack uf = cache

            uf.f = nlsolve_f(f, alg)
            uf.p = p
            uf.t = t

            J = jacobian(uf, uprev, integrator)
        end

        integrator.destats.njacs += 1

        if alg isa CompositeAlgorithm
            integrator.eigen_est = constvalue(opnorm(J, Inf))
        end
    end

    J
end

"""
    calc_J!(J, integrator, cache, next_step::Bool = false) -> J

Update the Jacobian object `J`.

If `integrator.f` has a custom Jacobian update function, then it will be called. Otherwise,
either automatic or finite differencing will be used depending on the `cache`.
If `next_step`, then it will evaluate the Jacobian at the next step.
"""
function calc_J!(J, integrator, cache, next_step::Bool = false)
    @unpack dt, t, uprev, f, p, alg = integrator
    if next_step
        t = t + dt
        uprev = integrator.u
    end

    if alg isa DAEAlgorithm
        if DiffEqBase.has_jac(f)
            duprev = integrator.duprev
            uf = cache.uf
            f.jac(J, duprev, uprev, p, uf.α * uf.invγdt, t)
        else
            @unpack du1, uf, jac_config = cache
            # using `dz` as temporary array
            x = cache.dz
            uf.t = t
            fill!(x, zero(eltype(x)))
            jacobian!(J, uf, x, du1, integrator, jac_config)
        end
    else
        if DiffEqBase.has_jac(f)
            f.jac(J, uprev, p, t)
        else
            @unpack du1, uf, jac_config = cache

            uf.f = nlsolve_f(f, alg)
            uf.t = t
            if !(p isa DiffEqBase.NullParameters)
                uf.p = p
            end

            jacobian!(J, uf, uprev, du1, integrator, jac_config)
        end
    end

    integrator.destats.njacs += 1

    if alg isa CompositeAlgorithm
        integrator.eigen_est = constvalue(opnorm(J, Inf))
    end

    return nothing
end

"""
    WOperator(mass_matrix,gamma,J[;transform=false])

A linear operator that represents the W matrix of an ODEProblem, defined as

```math
W = MM - \\gamma J
```

or, if `transform=true`:

```math
W = \\frac{1}{\\gamma}MM - J
```

where `MM` is the mass matrix (a regular `AbstractMatrix` or a `UniformScaling`),
`γ` is a real number proportional to the time step, and `J` is the Jacobian
operator (must be a `AbstractSciMLOperator`). A `WOperator` can also be
constructed using a `*DEFunction` directly as

    WOperator(f,gamma[;transform=false])

`f` needs to have a jacobian and `jac_prototype`, but the prototype does not need
to be a diffeq operator --- it will automatically be converted to one.

`WOperator` supports lazy `*` and `mul!` operations, the latter utilizing an
internal cache (can be specified in the constructor; default to regular `Vector`).
It supports all of `AbstractSciMLOperator`'s interface.
"""

"""
    islinearfunction(integrator) -> Tuple{Bool,Bool}

return the tuple `(is_linear_wrt_odealg, islinearodefunction)`.
"""
islinearfunction(integrator) = islinearfunction(integrator.f, integrator.alg)

"""
    islinearfunction(f, alg) -> Tuple{Bool,Bool}

return the tuple `(is_linear_wrt_odealg, islinearodefunction)`.
"""
function islinearfunction(f, alg)::Tuple{Bool, Bool}
    isode = f isa ODEFunction && islinear(f.f)
    islin = isode || (alg isa SplitAlgorithms && f isa SplitFunction && islinear(f.f1.f))
    return islin, isode
end

function do_newJW(integrator, alg, nlsolver, repeat_step)::NTuple{2, Bool}
    integrator.iter <= 1 && return true, true # at least one JW eval at the start
    repeat_step && return false, false
    islin, _ = islinearfunction(integrator)
    islin && return false, false # no further JW eval when it's linear
    !integrator.opts.adaptive && return true, true # Not adaptive will always refactorize
    alg isa DAEAlgorithm && return true, true
    isnewton(nlsolver) || return true, true
    isfirstcall(nlsolver) && return true, true
    isfs = isfirststage(nlsolver)
    iszero(nlsolver.fast_convergence_cutoff) && return isfs, isfs
    W_iγdt = inv(nlsolver.cache.W_γdt)
    iγdt = inv(nlsolver.γ * integrator.dt)
    smallstepchange = abs(iγdt / W_iγdt - 1) <= get_new_W_γdt_cutoff(nlsolver)
    jbad = nlsolver.status === TryAgain && smallstepchange
    errorfail = integrator.EEst > one(integrator.EEst)
    return jbad, (jbad || (!smallstepchange) || (isfs && errorfail))
end

@noinline _throwWJerror(W, J) = throw(DimensionMismatch("W: $(axes(W)), J: $(axes(J))"))
@noinline function _throwWMerror(W, mass_matrix)
    throw(DimensionMismatch("W: $(axes(W)), mass matrix: $(axes(mass_matrix))"))
end
@noinline function _throwJMerror(J, mass_matrix)
    throw(DimensionMismatch("J: $(axes(J)), mass matrix: $(axes(mass_matrix))"))
end

@noinline function calc_W(integrator, nlsolver, dtgamma, repeat_step, W_transform = false)
    @unpack t, uprev, p, f = integrator

    next_step = is_always_new(nlsolver)
    if next_step
        t = t + integrator.dt
        uprev = integrator.u
    end
    # Handle Rosenbrock has no nlsolver so passes cache directly
    cache = nlsolver isa OrdinaryDiffEqCache ? nlsolver : nlsolver.cache

    isdae = integrator.alg isa DAEAlgorithm
    if !isdae
        mass_matrix = integrator.f.mass_matrix
    end
    isarray = uprev isa AbstractArray
    # calculate W
    is_compos = integrator.alg isa CompositeAlgorithm
    islin, isode = islinearfunction(integrator)
    !isdae && update_coefficients!(mass_matrix, uprev, p, t)

    if islin
        J = isode ? f.f : f.f1.f # unwrap the Jacobian accordingly
        W = WOperator{false}(mass_matrix, dtgamma, J, uprev; transform = W_transform)
    elseif DiffEqBase.has_jac(f)
        J = f.jac(uprev, p, t)
        if typeof(J) <: StaticArray &&
           typeof(integrator.alg) <:
           Union{Rosenbrock23, Rodas4, Rodas4P, Rodas4P2, Rodas5, Rodas5P}
            W = W_transform ? J - mass_matrix * inv(dtgamma) :
                dtgamma * J - mass_matrix
        else
            if !isa(J, AbstractSciMLOperator) && (!isnewton(nlsolver) ||
                nlsolver.cache.W.J isa AbstractSciMLOperator)
                J = MatrixOperator(J)
            end
            W = WOperator{false}(mass_matrix, dtgamma, J, uprev, cache.W.jacvec;
                                 transform = W_transform)
        end
        integrator.destats.nw += 1
    else
        integrator.destats.nw += 1
        J = calc_J(integrator, cache, next_step)
        if isdae
            W = J
        else
            W_full = W_transform ? J - mass_matrix * inv(dtgamma) :
                     dtgamma * J - mass_matrix
            len = ArrayInterface.known_length(typeof(W_full))
            W = if W_full isa Number
                W_full
            elseif len !== nothing &&
                   typeof(integrator.alg) <:
                   Union{Rosenbrock23, Rodas4, Rodas4P, Rodas4P2, Rodas5, Rodas5P}
                StaticWOperator(W_full)
            else
                DiffEqBase.default_factorize(W_full)
            end
        end
    end
    (W isa WOperator && unwrap_alg(integrator, true) isa NewtonAlgorithm) &&
        (W = DiffEqBase.update_coefficients!(W, uprev, p, t)) # we will call `update_coefficients!` in NLNewton
    is_compos && (integrator.eigen_est = isarray ? constvalue(opnorm(J, Inf)) :
                            integrator.opts.internalnorm(J, t))
    return W
end

function calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step,
                                          W_transform)
    nlsolver = nothing
    # we need to skip calculating `W` when a step is repeated
    new_W = false
    if !repeat_step
        new_W = calc_W!(cache.W, integrator, nlsolver, cache, dtgamma, repeat_step,
                        W_transform)
    end
    calc_tderivative!(integrator, cache, dtd1, repeat_step)
    return new_W
end

# update W matrix (only used in Newton method)
function update_W!(integrator, cache, dtgamma, repeat_step, newJW = nothing)
    update_W!(cache.nlsolver, integrator, cache, dtgamma, repeat_step, newJW)
end

function update_W!(nlsolver::AbstractNLSolver,
                   integrator::SciMLBase.DEIntegrator{<:Any, true}, cache, dtgamma,
                   repeat_step::Bool, newJW = nothing)
    if isnewton(nlsolver)
        calc_W!(get_W(nlsolver), integrator, nlsolver, cache, dtgamma, repeat_step, true,
                newJW)
    end
    nothing
end

function update_W!(nlsolver::AbstractNLSolver,
                   integrator::SciMLBase.DEIntegrator{<:Any, false}, cache, dtgamma,
                   repeat_step::Bool, newJW = nothing)
    if isnewton(nlsolver)
        isdae = integrator.alg isa DAEAlgorithm
        new_jac, new_W = true, true
        if isdae && new_jac
            lcache = nlsolver.cache
            lcache.uf.α = nlsolver.α
            lcache.uf.invγdt = inv(dtgamma)
            lcache.uf.tmp = @. nlsolver.tmp
            lcache.uf.uprev = @. integrator.uprev
        end
        nlsolver.cache.W = calc_W(integrator, nlsolver, dtgamma, repeat_step, true)
        #TODO: jacobian reuse for oop
        new_jac && (nlsolver.cache.J_t = integrator.t)
        set_new_W!(nlsolver, new_W)
        if new_jac && isdae
            set_W_γdt!(nlsolver, nlsolver.α * inv(dtgamma))
        elseif new_W && !isdae
            set_W_γdt!(nlsolver, dtgamma)
        end
    end
    nothing
end

SciMLBase.isinplace(::WOperator{IIP}, i) where {IIP} = IIP
set_gamma!(W::WOperator, gamma) = (W.gamma = gamma; W)

function SciMLOperators.update_coefficients!(J::SparseDiffTools.JacVec,u,p,t)
  copyto!(J.x,u)
  J.f.t = t
  J.f.p = p
end

function build_J_W(alg, u, p, t, dt, f, γ, IIP, transform=true)

    # special consiederations for splitfunction, dynamicalodefunction?

    ## FLAGS
    # IIP
    # concrete_jac(alg)
    # LinearSolve.needs_concrete_A(linsolve)

    γ = if γ isa ScalarOperator
        γ
    else
        update_func = DEFAULT_UPDATE_FUNC # TODO
        ScalarOperator(γ; update_func = update_func)
    end

    mm = f.mass_matrix
    M = if mm isa AbstractSciMLOperator
        mm
    elseif mm isa AbstractMatrix
        MatrixOperator(mm)
    end

    J = if !isa(f.jac, Nothing)
        if f.jac isa AbstractSciMLOperator
            f.jac
        elseif f.jac isa AbstractMatrix
            MatrixOperator(f.jac)
        elseif f.jac isa AbstractSparseMatrix
            # TODO


            # Workaround https://github.com/JuliaSparse/SparseArrays.jl/issues/190
            # Hopefully `rand()` does not match any value in the array (prob ~ 0, with a check)
            # Then `one` is required since gamma is zero
            # Otherwise this will not pick up the union sparsity pattern
            # But instead drop the runtime zeros (i.e. all values) of the AJ pattern!
            Jn = nonzeros(f.jac)
            x = rand()
            @assert all(!isequal(x), Jn)
            fill!(AJn, rand())
            !(_concrete_form, false) # safety measure, throw singular error if not filled
        else
            f.jac
        end
    elseif !isa(f.jac_prototype, Nothing)
        f.jac_prototype |> deepcopy |> MatrixOperator
    elseif !isa(f.jvp, Nothing)
        jvp isa FunctionOperator ? f.jvp : FunctionOperator(f.jvp, u, u, p=p, t=t)
    else
        islin, isode = islinearfunction(f, alg)
        _f = islin ? (isode ? f.f : f.f1.f) : f
        SparseDiffTools.JacVec(UJacobianWrapper(_f, t, p),
                               copy(u),
                               OrdinaryDiffEqTag(),
                               autodiff = alg_autodiff(alg))
    end

    M = cache_operator(M, u)
    J = cache_operator(J, u)

    W = if transform
        (1 / γ) * M - J
    else
        M - γ * J 
    end

    W = cache_operator(W, u)

    # update_coefficients!(W, u, p, t)

    J, W
end

function build_J_W(alg, u, uprev, p, t, dt, f::F, ::Type{uEltypeNoUnits},
                   ::Val{IIP}) where {IIP, uEltypeNoUnits, F}
    islin, isode = islinearfunction(f, alg)
    if f.jac_prototype isa AbstractSciMLOperator
        W = WOperator{IIP}(f, u, dt)
        J = W.J
    elseif IIP && f.jac_prototype !== nothing && concrete_jac(alg) === nothing &&
           (alg.linsolve === nothing ||
            alg.linsolve !== nothing &&
            LinearSolve.needs_concrete_A(alg.linsolve))

        # If factorization, then just use the jac_prototype
        J = similar(f.jac_prototype)
        W = similar(J)
    elseif (IIP && (concrete_jac(alg) === nothing || !concrete_jac(alg)) &&
            alg.linsolve !== nothing &&
            !LinearSolve.needs_concrete_A(alg.linsolve))
        # If the user has chosen GMRES but no sparse Jacobian, assume that the dense
        # Jacobian is a bad idea and create a fully matrix-free solver. This can
        # be overriden with concrete_jac.

        _f = islin ? (isode ? f.f : f.f1.f) : f
        jacvec = SparseDiffTools.JacVec(UJacobianWrapper(_f, t, p), copy(u),
                                        OrdinaryDiffEqTag(), autodiff = alg_autodiff(alg))
        J = jacvec
        W = WOperator{IIP}(f.mass_matrix, dt, J, u, jacvec)

    elseif alg.linsolve !== nothing && !LinearSolve.needs_concrete_A(alg.linsolve) ||
           concrete_jac(alg) !== nothing && concrete_jac(alg)
        # The linear solver does not need a concrete Jacobian, but the user has
        # asked for one. This will happen when the Jacobian is used in the preconditioner
        # Thus setup JacVec and a concrete J, using sparsity when possible
        _f = islin ? (isode ? f.f : f.f1.f) : f
        J = if f.jac_prototype === nothing
            ArrayInterface.zeromatrix(u)
        else
            deepcopy(f.jac_prototype)
        end
        jacvec = SparseDiffTools.JacVec(UJacobianWrapper(_f, t, p), copy(u),
                                        OrdinaryDiffEqTag(), autodiff = alg_autodiff(alg))
        W = WOperator{IIP}(f.mass_matrix, dt, J, u, jacvec)

    elseif islin || (!IIP && DiffEqBase.has_jac(f))
        J = islin ? (isode ? f.f : f.f1.f) : f.jac(uprev, p, t) # unwrap the Jacobian accordingly
        if !isa(J, AbstractSciMLOperator)
            J = MatrixOperator(J)
        end
        W = WOperator{IIP}(f.mass_matrix, dt, J, u)
    else
        J = if f.jac_prototype === nothing
            ArrayInterface.zeromatrix(u)
        else
            deepcopy(f.jac_prototype)
        end
        isdae = alg isa DAEAlgorithm
        W = if isdae
            J
        elseif IIP
            similar(J)
        else
            len = ArrayInterface.known_length(typeof(J))
            if len !== nothing &&
                typeof(alg) <:
                Union{Rosenbrock23, Rodas4, Rodas4P, Rodas4P2, Rodas5, Rodas5P}
                StaticWOperator(J, false)
            else
                ArrayInterface.lu_instance(J)
            end
        end
    end
    return J, W
end

build_uf(alg, nf, t, p, ::Val{true}) = UJacobianWrapper(nf, t, p)
build_uf(alg, nf, t, p, ::Val{false}) = UDerivativeWrapper(nf, t, p)
