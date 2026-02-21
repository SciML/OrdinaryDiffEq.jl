using SciMLOperators: StaticWOperator, WOperator

function calc_tderivative!(integrator, cache, dtd1, repeat_step)
    return @inbounds begin
        (; t, dt, uprev, u, f, p) = integrator
        (; du2, fsalfirst, dT, tf, linsolve_tmp) = cache

        # Time derivative
        if !repeat_step # skip calculation if step is repeated
            if SciMLBase.has_tgrad(f)
                f.tgrad(dT, uprev, p, t)
            else
                tf.uprev = uprev
                tf.p = p
                alg = unwrap_alg(integrator, true)

                autodiff_alg = ADTypes.dense_ad(gpu_safe_autodiff(alg_autodiff(alg), u))

                # Convert t to eltype(dT) if using ForwardDiff, to make FunctionWrappers work
                t = autodiff_alg isa AutoForwardDiff ? convert(eltype(dT), t) : t

                grad_config_tup = cache.grad_config

                if autodiff_alg isa AutoFiniteDiff
                    grad_config = diffdir(integrator) > 0 ? grad_config_tup[1] :
                        grad_config_tup[2]
                else
                    grad_config = grad_config_tup[1]
                end

                if integrator.iter == 1
                    try
                        DI.derivative!(
                            tf, linsolve_tmp, dT, grad_config, autodiff_alg, t
                        )
                    catch e
                        throw(FirstAutodiffTgradError(e))
                    end
                else
                    DI.derivative!(tf, linsolve_tmp, dT, grad_config, autodiff_alg, t)
                end

                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            end
        end

        @.. broadcast = false linsolve_tmp = fsalfirst + dtd1 * dT
    end
end

function calc_tderivative(integrator, cache)
    (; t, dt, uprev, u, f, p, alg) = integrator

    # Time derivative
    if SciMLBase.has_tgrad(f)
        dT = f.tgrad(uprev, p, t)
    else
        tf = cache.tf
        tf.u = uprev
        tf.p = p

        autodiff_alg = ADTypes.dense_ad(gpu_safe_autodiff(alg_autodiff(alg), u))

        if alg_autodiff isa AutoFiniteDiff
            autodiff_alg = SciMLBase.@set autodiff_alg.dir = diffdir(integrator)
        end

        if integrator.iter == 1
            try
                dT = DI.derivative(tf, autodiff_alg, t)
            catch e
                throw(FirstAutodiffTgradError(e))
            end
        else
            dT = DI.derivative(tf, autodiff_alg, t)
        end

        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
    return dT
end

"""
    calc_J(integrator, cache, next_step::Bool = false)

Return a new Jacobian object.

If `integrator.f` has a custom Jacobian update function, then it will be called. Otherwise,
either automatic or finite differencing will be used depending on the `uf` object of the
cache. If `next_step`, then it will evaluate the Jacobian at the next step.
"""
function calc_J(integrator, cache, next_step::Bool = false)
    (; dt, t, uprev, f, p, alg) = integrator
    if next_step
        t = t + dt
        uprev = integrator.u
    end

    method = if SciMLBase.has_jac(f)
        "user-provided"
    else
        if hasproperty(cache, :jac_config) && cache.jac_config !== nothing
            "autodiff"
        else
            "finite-diff"
        end
    end

    @SciMLMessage(
        lazy"Computing Jacobian at t = $(t) using $(method)",
        integrator.opts.verbose, :jacobian_update
    )

    if alg isa DAEAlgorithm
        if SciMLBase.has_jac(f)
            duprev = integrator.duprev
            uf = cache.uf
            J = f.jac(duprev, uprev, p, uf.α * uf.invγdt, t)
        else
            (; uf) = cache
            x = zero(uprev)
            J = jacobian(uf, x, integrator)
        end
    else
        if SciMLBase.has_jac(f)
            J = f.jac(uprev, p, t)
        else
            (; uf) = cache

            uf.f = nlsolve_f(f, alg)
            uf.p = p
            uf.t = t
            J = jacobian(uf, uprev, integrator)
        end

        if alg isa CompositeAlgorithm
            integrator.eigen_est = constvalue(opnorm(J, Inf))
        end
    end

    integrator.stats.njacs += 1
    return J
end

"""
    calc_J!(J, integrator, cache, next_step::Bool = false) -> J

Update the Jacobian object `J`.

If `integrator.f` has a custom Jacobian update function, then it will be called. Otherwise,
either automatic or finite differencing will be used depending on the `cache`.
If `next_step`, then it will evaluate the Jacobian at the next step.
"""
function calc_J!(J, integrator, cache, next_step::Bool = false)
    (; dt, t, uprev, f, p, alg) = integrator
    if next_step
        t = t + dt
        uprev = integrator.u
    end

    if alg isa DAEAlgorithm
        if SciMLBase.has_jac(f)
            duprev = integrator.duprev
            uf = cache.uf
            # need to do some jank here to account for sparsity pattern of W
            # https://github.com/SciML/OrdinaryDiffEq.jl/issues/2653

            # we need to set all nzval to a non-zero number
            # otherwise in the following line any zero gets interpreted as a structural zero
            if !isnothing(integrator.f.jac_prototype) &&
                    is_sparse_csc(integrator.f.jac_prototype)
                set_all_nzval!(integrator.f.jac_prototype, true)
                J .= true .* integrator.f.jac_prototype
                set_all_nzval!(J, false)
                f.jac(J, duprev, uprev, p, uf.α * uf.invγdt, t)
            else
                f.jac(J, duprev, uprev, p, uf.α * uf.invγdt, t)
            end
        else
            (; du1, uf, jac_config) = cache
            # using `dz` as temporary array
            x = cache.dz
            uf.t = t
            fill!(x, zero(eltype(x)))
            jacobian!(J, uf, x, du1, integrator, jac_config)
        end
    else
        if SciMLBase.has_jac(f)
            # need to do some jank here to account for sparsity pattern of W
            # https://github.com/SciML/OrdinaryDiffEq.jl/issues/2653

            # we need to set all nzval to a non-zero number
            # otherwise in the following line any zero gets interpreted as a structural zero
            if !isnothing(integrator.f.jac_prototype) &&
                    is_sparse_csc(integrator.f.jac_prototype)
                set_all_nzval!(integrator.f.jac_prototype, true)
                J .= true .* integrator.f.jac_prototype
                set_all_nzval!(J, false)
                f.jac(J, uprev, p, t)
            else
                f.jac(J, uprev, p, t)
            end
        else
            (; du1, uf, jac_config) = cache
            uf.f = nlsolve_f(f, alg)
            uf.t = t
            if !(p isa SciMLBase.NullParameters)
                uf.p = p
            end
            jacobian!(J, uf, uprev, du1, integrator, jac_config)
        end
    end

    if alg isa CompositeAlgorithm
        integrator.eigen_est = constvalue(opnorm(J, Inf))
    end

    integrator.stats.njacs += 1
    return nothing
end

"""
    islinearfunction(integrator) -> Tuple{Bool,Bool}

return the tuple `(is_linear_wrt_odealg, islinearodefunction)`.
"""
islinearfunction(integrator) = islinearfunction(integrator.f, integrator.alg)

"""
    islinearfunction(f, alg) -> Tuple{Bool,Bool}

return the tuple `(is_linear_wrt_odealg, islinearodefunction)`.
"""
function islinearfunction(f::F, alg)::Tuple{Bool, Bool} where {F}
    isode = f isa ODEFunction && islinear(f.f)
    islin = isode || (issplit(alg) && f isa SplitFunction && islinear(f.f1.f))
    return islin, isode
end

function do_newJW(integrator, alg, nlsolver, repeat_step)::NTuple{2, Bool}
    integrator.iter <= 1 && return true, true # at least one JW eval at the start
    repeat_step && return false, false
    islin, _ = islinearfunction(integrator)
    islin && return false, false # no further JW eval when it's linear
    !integrator.opts.adaptive && return true, true # Not adaptive will always refactorize
    errorfail = integrator.EEst > one(integrator.EEst)
    # TODO: add `isJcurrent` support for Rosenbrock solvers
    if !isnewton(nlsolver)
        isfreshJ = !(integrator.alg isa CompositeAlgorithm) &&
            (integrator.iter > 1 && errorfail && !integrator.u_modified)
        return !isfreshJ, true
    end
    isfirstcall(nlsolver) && return true, true
    isfs = isfirststage(nlsolver)
    isfreshJ = isJcurrent(nlsolver, integrator) && !integrator.u_modified
    iszero(nlsolver.fast_convergence_cutoff) && return isfs && !isfreshJ, isfs
    isdae = alg isa DAEAlgorithm
    if !isdae
        mm = integrator.f.mass_matrix
        is_varying_mm = !isconstant(mm)
    end
    if isfreshJ
        jbad = false
        smallstepchange = true
    else
        if isdae
            # IDA-style cj ratio test for DAE solvers.
            # For DAE, W_γdt stores α / (γ * dt) from the last Jacobian eval.
            # Compare current cj = α / (γ * dt) with stored value.
            current_cj = nlsolver.α * inv(nlsolver.γ * integrator.dt)
            old_cj = nlsolver.cache.W_γdt
            smallstepchange = abs(current_cj / old_cj - 1) <=
                get_new_W_γdt_cutoff(nlsolver)
        else
            W_iγdt = inv(nlsolver.cache.W_γdt)
            iγdt = inv(nlsolver.γ * integrator.dt)
            smallstepchange = abs(iγdt / W_iγdt - 1) <=
                get_new_W_γdt_cutoff(nlsolver)
        end
        jbad = nlsolver.status === TryAgain && smallstepchange
    end
    wbad = (!smallstepchange) || (isfs && errorfail) || nlsolver.status === Divergence
    if isdae
        # For DAE solvers, W = J (no separate I - γ*dt*J transformation).
        # Recomputing W without J is a no-op, so new_jac and new_W must match.
        need_recompute = jbad || wbad
        return need_recompute, need_recompute
    else
        return jbad, (is_varying_mm || jbad || wbad)
    end
end

@noinline _throwWJerror(W, J) = throw(DimensionMismatch("W: $(axes(W)), J: $(axes(J))"))
@noinline function _throwWMerror(W, mass_matrix)
    throw(DimensionMismatch("W: $(axes(W)), mass matrix: $(axes(mass_matrix))"))
end
@noinline function _throwJMerror(J, mass_matrix)
    throw(DimensionMismatch("J: $(axes(J)), mass matrix: $(axes(mass_matrix))"))
end

function jacobian2W!(
        W::AbstractMatrix, mass_matrix, dtgamma::Number, J::AbstractMatrix
    )::Nothing
    # check size and dimension
    iijj = axes(W)
    @boundscheck (iijj == axes(J) && length(iijj) == 2) || _throwWJerror(W, J)
    mass_matrix isa UniformScaling ||
        @boundscheck axes(mass_matrix) == axes(W) || _throwWMerror(W, mass_matrix)
    @inbounds begin
        invdtgamma = inv(dtgamma)
        if mass_matrix isa UniformScaling
            copyto!(W, J)
            idxs = diagind(W)
            λ = -mass_matrix.λ
            if ArrayInterface.fast_scalar_indexing(J) &&
                    ArrayInterface.fast_scalar_indexing(W)
                @inbounds for i in 1:size(J, 1)
                    W[i, i] = muladd(λ, invdtgamma, J[i, i])
                end
            else
                @.. broadcast = false @view(W[idxs]) = muladd(λ, invdtgamma, @view(J[idxs]))
            end
        else
            @.. broadcast = false W = muladd(-mass_matrix, invdtgamma, J)
        end
    end
    return nothing
end

function jacobian2W!(W::Matrix, mass_matrix, dtgamma::Number, J::Matrix)::Nothing
    # check size and dimension
    iijj = axes(W)
    @boundscheck (iijj == axes(J) && length(iijj) == 2) || _throwWJerror(W, J)
    mass_matrix isa UniformScaling ||
        @boundscheck axes(mass_matrix) == axes(W) || _throwWMerror(W, mass_matrix)
    @inbounds begin
        invdtgamma = inv(dtgamma)
        if mass_matrix isa UniformScaling
            copyto!(W, J)
            idxs = diagind(W)
            λ = -mass_matrix.λ
            @inbounds for i in 1:size(J, 1)
                W[i, i] = muladd(λ, invdtgamma, J[i, i])
            end
        else
            @inbounds @simd ivdep for i in eachindex(W)
                W[i] = muladd(-mass_matrix[i], invdtgamma, J[i])
            end
        end
    end
    return nothing
end

function jacobian2W(mass_matrix, dtgamma::Number, J::AbstractMatrix)
    # check size and dimension
    mass_matrix isa UniformScaling ||
        @boundscheck axes(mass_matrix) == axes(J) || _throwJMerror(J, mass_matrix)
    @inbounds begin
        invdtgamma = inv(dtgamma)
        if mass_matrix isa UniformScaling
            λ = -mass_matrix.λ
            W = J + (λ * invdtgamma) * I
        else
            W = muladd(-mass_matrix, invdtgamma, J)
        end
    end
    return W
end

is_always_new(alg) = isdefined(alg, :always_new) ? alg.always_new : false

function calc_W!(
        W, integrator, nlsolver::Union{Nothing, AbstractNLSolver}, cache, dtgamma,
        repeat_step, newJW = nothing
    )
    (; t, dt, uprev, u, f, p) = integrator
    lcache = nlsolver === nothing ? cache : nlsolver.cache
    next_step = is_always_new(nlsolver)
    if next_step
        t = t + integrator.dt
        uprev = integrator.u
    end

    (; J) = lcache
    isdae = integrator.alg isa DAEAlgorithm
    alg = unwrap_alg(integrator, true)
    mass_matrix = nothing
    if !isdae
        mass_matrix = integrator.f.mass_matrix
    end
    is_compos = integrator.alg isa CompositeAlgorithm

    # handle Wfact
    if SciMLBase.has_Wfact_t(f)
        f.Wfact_t(W, u, p, dtgamma, t)
        isnewton(nlsolver) && set_W_γdt!(nlsolver, dtgamma)
        is_compos && (
            integrator.eigen_est = constvalue(opnorm(LowerTriangular(W), Inf)) +
                inv(dtgamma)
        ) # TODO: better estimate
        # It's equivalent with evaluating a new Jacobian, but not a new W,
        # because we won't call `lu!`, and the iteration matrix is fresh.
        return (true, false)
    end

    # check if we need to update J or W
    if newJW === nothing
        new_jac, new_W = do_newJW(integrator, alg, nlsolver, repeat_step)
    else
        new_jac, new_W = newJW
    end

    if new_jac && isnewton(lcache)
        lcache.J_t = t
        if isdae
            lcache.uf.α = nlsolver.α
            lcache.uf.invγdt = inv(dtgamma)
            lcache.uf.tmp = nlsolver.tmp
        end
    end

    # calculate W
    if W isa WOperator
        if isnewton(nlsolver)
            # we will call `update_coefficients!` for u/p/t in NLNewton
            update_coefficients!(W; gamma = dtgamma)
        else
            update_coefficients!(W, uprev, p, t; gamma = dtgamma)
        end
        if W.J !== nothing && !(W.J isa AbstractSciMLOperator)
            islin, isode = islinearfunction(integrator)
            islin ? (J = isode ? f.f : f.f1.f) :
                (new_jac && (calc_J!(W.J, integrator, lcache, next_step)))
            new_W && !isdae &&
                jacobian2W!(W._concrete_form, mass_matrix, dtgamma, J)
        end
    elseif W isa AbstractSciMLOperator && !(W isa StaticWOperator)
        update_coefficients!(W, uprev, p, t; gamma = dtgamma)
    else # concrete W using jacobian from `calc_J!`
        islin, isode = islinearfunction(integrator)
        islin ? (J = isode ? f.f : f.f1.f) :
            (new_jac && (calc_J!(J, integrator, lcache, next_step)))
        new_W && !isdae && jacobian2W!(W, mass_matrix, dtgamma, J)
    end
    if isnewton(nlsolver)
        set_new_W!(nlsolver, new_W)
        if new_jac && isdae
            set_W_γdt!(nlsolver, nlsolver.α * inv(dtgamma))
        elseif new_W && !isdae
            set_W_γdt!(nlsolver, dtgamma)
        end
    end

    new_W && (integrator.stats.nw += 1)
    return new_jac, new_W
end

@noinline function calc_W(integrator, nlsolver, dtgamma, repeat_step)
    (; t, uprev, p, f) = integrator

    next_step = is_always_new(nlsolver)
    if next_step
        t = t + integrator.dt
        uprev = integrator.u
    end
    # Handle Rosenbrock has no nlsolver so passes cache directly
    cache = nlsolver isa OrdinaryDiffEqCache ? nlsolver : nlsolver.cache

    isdae = integrator.alg isa DAEAlgorithm
    mass_matrix = nothing
    if !isdae
        mass_matrix = integrator.f.mass_matrix
    end
    isarray = uprev isa AbstractArray
    # calculate W
    is_compos = integrator.alg isa CompositeAlgorithm
    islin, isode = islinearfunction(integrator)
    !isdae && update_coefficients!(mass_matrix, uprev, p, t)

    J = nothing
    if cache.W isa StaticWOperator
        integrator.stats.nw += 1
        J = calc_J(integrator, cache, next_step)
        W = StaticWOperator(J - mass_matrix * inv(dtgamma))
    elseif cache.W isa WOperator
        integrator.stats.nw += 1
        J = if islin
            isode ? f.f : f.f1.f
        else
            calc_J(integrator, cache, next_step)
        end
        W = WOperator{false}(mass_matrix, dtgamma, J, uprev, cache.W.jacvec)
    elseif cache.W isa AbstractSciMLOperator
        W = update_coefficients(cache.W, uprev, p, t; gamma = dtgamma)
    else
        integrator.stats.nw += 1
        J = islin ? isode ? f.f : f.f1.f : calc_J(integrator, cache, next_step)
        if isdae
            W = J
        else
            W = J - mass_matrix * inv(dtgamma)

            if !isa(W, Number)
                W = DiffEqBase.default_factorize(W)
            end
        end
    end
    is_compos && (
        integrator.eigen_est = isarray ? constvalue(opnorm(J, Inf)) :
            integrator.opts.internalnorm(J, t)
    )
    return W
end

function calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step)
    nlsolver = nothing
    # we need to skip calculating `J` and `W` when a step is repeated
    new_jac = new_W = false
    if !repeat_step
        new_jac, new_W = calc_W!(
            cache.W, integrator, nlsolver, cache, dtgamma, repeat_step
        )
    end
    # If the Jacobian is not updated, we won't have to update ∂/∂t either.
    calc_tderivative!(integrator, cache, dtd1, repeat_step || !new_jac)
    return new_W
end

# update W matrix (only used in Newton method)
function update_W!(integrator, cache, dtgamma, repeat_step, newJW = nothing)
    return update_W!(cache.nlsolver, integrator, cache, dtgamma, repeat_step, newJW)
end

function update_W!(
        nlsolver::AbstractNLSolver,
        integrator::SciMLBase.DEIntegrator{<:Any, true}, cache, dtgamma,
        repeat_step::Bool, newJW = nothing
    )
    if isnewton(nlsolver)
        new_jac, new_W = calc_W!(
            get_W(nlsolver), integrator, nlsolver, cache, dtgamma, repeat_step,
            newJW
        )
        if new_W
            @SciMLMessage(
                lazy"W matrix factorized: dtgamma = $(dtgamma), new_jac = $(new_jac)",
                integrator.opts.verbose, :w_factorization
            )
        end
    end
    return nothing
end

function update_W!(
        nlsolver::AbstractNLSolver,
        integrator::SciMLBase.DEIntegrator{<:Any, false}, cache, dtgamma,
        repeat_step::Bool, newJW = nothing
    )
    if isnewton(nlsolver)
        isdae = integrator.alg isa DAEAlgorithm
        if newJW === nothing
            new_jac, new_W = do_newJW(integrator, integrator.alg, nlsolver, repeat_step)
        else
            new_jac, new_W = newJW
        end
        if isdae && new_jac
            lcache = nlsolver.cache
            lcache.uf.α = nlsolver.α
            lcache.uf.invγdt = inv(dtgamma)
            lcache.uf.tmp = @. nlsolver.tmp
            lcache.uf.uprev = @. integrator.uprev
        end
        if new_W
            nlsolver.cache.W = calc_W(integrator, nlsolver, dtgamma, repeat_step)
        end
        new_jac && (nlsolver.cache.J_t = integrator.t)
        set_new_W!(nlsolver, new_W)
        if new_jac && isdae
            set_W_γdt!(nlsolver, nlsolver.α * inv(dtgamma))
        elseif new_W && !isdae
            set_W_γdt!(nlsolver, dtgamma)
        end
        if new_W
            @SciMLMessage(
                lazy"W matrix factorized: dtgamma = $(dtgamma), new_jac = $(new_jac)",
                integrator.opts.verbose, :w_factorization
            )
        end
    end
    return nothing
end

function build_J_W(
        alg, u, uprev, p, t, dt, f::F, jac_config, ::Type{uEltypeNoUnits},
        ::Val{IIP}
    ) where {IIP, uEltypeNoUnits, F}
    # TODO - make J, W AbstractSciMLOperators (lazily defined with scimlops functionality)
    # TODO - if jvp given, make it SciMLOperators.FunctionOperator
    # TODO - make mass matrix a SciMLOperator so it can be updated with time. Default to IdentityOperator
    islin, isode = islinearfunction(f, alg)
    if isdefined(f, :W_prototype) && (f.W_prototype isa AbstractSciMLOperator)
        # We use W_prototype when it is provided as a SciMLOperator, and in this case we require jac_prototype to be a SciMLOperator too.
        if !(f.jac_prototype isa AbstractSciMLOperator)
            error("SciMLOperator for W_prototype only supported when jac_prototype is a SciMLOperator, but got $(typeof(f.jac_prototype))")
        end
        W = f.W_prototype
        J = f.jac_prototype
    elseif f.jac_prototype isa AbstractSciMLOperator
        J = deepcopy(f.jac_prototype)
        if J isa AbstractMatrix
            @assert SciMLBase.has_jac(f) "f needs to have an associated jacobian"
            J = MatrixOperator(J; update_func! = f.jac)
        end
        W = WOperator{IIP}(f.mass_matrix, promote(t, dt)[2], J, _vec(u))
    elseif islin
        J = isode ? f.f : f.f1.f # unwrap the Jacobian accordingly
        W = WOperator{IIP}(f.mass_matrix, dt, J, _vec(u))
    elseif IIP && f.jac_prototype !== nothing && concrete_jac(alg) === nothing &&
            (alg.linsolve === nothing || LinearSolve.needs_concrete_A(alg.linsolve))

        # If factorization, then just use the jac_prototype
        J = similar(f.jac_prototype)
        W = similar(J)
    elseif (
            IIP && (concrete_jac(alg) === nothing || !concrete_jac(alg)) &&
                alg.linsolve !== nothing &&
                !LinearSolve.needs_concrete_A(alg.linsolve)
        )
        # If the user has chosen GMRES but no sparse Jacobian, assume that the dense
        # Jacobian is a bad idea and create a fully matrix-free solver. This can
        # be overridden with concrete_jac.
        jacvec = JVPCache(f, copy(u), u, p, t, autodiff = alg_autodiff(alg))

        J = jacvec
        W = WOperator{IIP}(f.mass_matrix, promote(t, dt)[2], J, _vec(u), jacvec)
    elseif alg.linsolve !== nothing && !LinearSolve.needs_concrete_A(alg.linsolve) ||
            concrete_jac(alg) !== nothing && concrete_jac(alg)
        # The linear solver does not need a concrete Jacobian, but the user has
        # asked for one. This will happen when the Jacobian is used in the preconditioner
        # Thus setup JacVec and a concrete J, using sparsity when possible
        _f = islin ? (isode ? f.f : f.f1.f) : f
        J = if f.jac_prototype === nothing
            if alg_autodiff(alg) isa AutoSparse
                if isnothing(f.sparsity)
                    !isnothing(jac_config) ?
                        convert.(
                            eltype(u), sparsity_pattern(jac_config[1])
                        ) :
                        spzeros(eltype(u), length(u), length(u))
                elseif eltype(f.sparsity) == Bool
                    convert.(eltype(u), f.sparsity)
                else
                    f.sparsity
                end
            else
                ArrayInterface.zeromatrix(u)
            end
        else
            deepcopy(f.jac_prototype)
        end
        W = if J isa StaticMatrix
            StaticWOperator(J, false)
        else
            jacvec = JVPCache(f, copy(u), u, p, t, autodiff = alg_autodiff(alg))

            WOperator{IIP}(f.mass_matrix, promote(t, dt)[2], J, _vec(u), jacvec)
        end
    else
        J = if !IIP && SciMLBase.has_jac(f)
            if f isa DAEFunction
                f.jac(uprev, uprev, p, one(t), t)
            else
                f.jac(uprev, p, t)
            end
        elseif f.jac_prototype === nothing
            if alg_autodiff(alg) isa AutoSparse
                if isnothing(f.sparsity)
                    !isnothing(jac_config) ? convert.(eltype(u), sparsity_pattern(jac_config[1])) :
                        spzeros(eltype(u), length(u), length(u))
                elseif eltype(f.sparsity) == Bool
                    convert.(eltype(u), f.sparsity)
                else
                    f.sparsity
                end
            else
                ArrayInterface.zeromatrix(u)
            end
        else
            deepcopy(f.jac_prototype)
        end
        W = if alg isa DAEAlgorithm
            J
        elseif IIP
            similar(J)
        elseif J isa StaticMatrix
            StaticWOperator(J, false)
        else
            ArrayInterface.lu_instance(J)
        end
    end
    return J, W
end

build_uf(alg, nf, t, p, ::Val{true}) = UJacobianWrapper(nf, t, p)
build_uf(alg, nf, t, p, ::Val{false}) = UDerivativeWrapper(nf, t, p)

function LinearSolve.init_cacheval(
        alg::LinearSolve.DefaultLinearSolver, A::WOperator, b, u,
        Pl, Pr,
        maxiters::Int, abstol, reltol, verbose::LinearVerbosity,
        assumptions::OperatorAssumptions
    )
    return LinearSolve.init_cacheval(
        alg, A.J, b, u, Pl, Pr,
        maxiters::Int, abstol, reltol, verbose::LinearVerbosity,
        assumptions::OperatorAssumptions
    )
end

for alg in [
        LinearSolve.AppleAccelerateLUFactorization,
        LinearSolve.BunchKaufmanFactorization,
        LinearSolve.CHOLMODFactorization,
        LinearSolve.CholeskyFactorization,
        LinearSolve.CudaOffloadFactorization,
        LinearSolve.DiagonalFactorization,
        LinearSolve.FastLUFactorization,
        LinearSolve.FastQRFactorization,
        LinearSolve.GenericFactorization,
        LinearSolve.GenericLUFactorization,
        LinearSolve.KLUFactorization,
        LinearSolve.LDLtFactorization,
        LinearSolve.LUFactorization,
        LinearSolve.MKLLUFactorization,
        LinearSolve.MetalLUFactorization,
        LinearSolve.NormalBunchKaufmanFactorization,
        LinearSolve.NormalCholeskyFactorization,
        LinearSolve.QRFactorization,
        LinearSolve.RFLUFactorization,
        LinearSolve.SVDFactorization,
        LinearSolve.SimpleLUFactorization,
        LinearSolve.SparspakFactorization,
        LinearSolve.UMFPACKFactorization,
    ]
    @eval function LinearSolve.init_cacheval(
            alg::$alg, A::WOperator, b, u, Pl, Pr,
            maxiters::Int, abstol, reltol, verbose::LinearVerbosity,
            assumptions::OperatorAssumptions
        )
        return LinearSolve.init_cacheval(
            alg, A.J, b, u, Pl, Pr,
            maxiters::Int, abstol, reltol, verbose::LinearVerbosity,
            assumptions::OperatorAssumptions
        )
    end
end

function resize_J_W!(cache, integrator, i)
    (isdefined(cache, :J) && isdefined(cache, :W)) || return

    (; f) = integrator

    if cache.W isa WOperator
        nf = nlsolve_f(f, integrator.alg)
        islin = f isa Union{ODEFunction, SplitFunction} && islinear(nf.f)
        if !islin
            if cache.J isa AbstractSciMLOperator
                resize_JVPCache!(
                    cache.J, f, cache.du1, integrator.u, alg_autodiff(integrator.alg)
                )
            elseif f.jac_prototype !== nothing
                J = similar(f.jac_prototype, i, i)
                J = MatrixOperator(J; update_func! = f.jac)
            end
            if cache.W.jacvec isa AbstractSciMLOperator
                resize_JVPCache!(
                    cache.W.jacvec, f, cache.du1, integrator.u,
                    alg_autodiff(integrator.alg)
                )
            end
            cache.W = WOperator{SciMLBase.isinplace(integrator.sol.prob)}(
                f.mass_matrix,
                integrator.dt,
                cache.J,
                integrator.u,
                cache.W.jacvec
            )
            cache.J = cache.W.J
        end
    else
        if cache.J !== nothing
            cache.J = similar(cache.J, i, i)
        end
        cache.W = similar(cache.W, i, i)
    end

    return nothing
end

getsize(::Val{N}) where {N} = N
getsize(N::Integer) = N
