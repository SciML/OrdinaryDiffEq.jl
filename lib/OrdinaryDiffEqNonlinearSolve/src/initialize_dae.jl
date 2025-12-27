# Efficient algebraic variable/equation detection for sparse mass matrices.
# O(nnz) instead of O(n²) for sparse matrices.
"""
    find_algebraic_vars_eqs(M::SparseMatrixCSC)

Find algebraic variables (zero columns) and algebraic equations (zero rows) from mass matrix.
Returns (algebraic_vars::Vector{Bool}, algebraic_eqs::Vector{Bool}).

For sparse matrices, uses O(nnz) traversal of CSC structure instead of O(n²) iteration.
"""
function find_algebraic_vars_eqs(M::SparseMatrixCSC)
    n_cols = size(M, 2)
    n_rows = size(M, 1)

    # Initialize all as algebraic (true = zero column/row)
    algebraic_vars = fill(true, n_cols)
    algebraic_eqs = fill(true, n_rows)

    # Mark columns/rows with non-zero values as differential (false)
    @inbounds for j in 1:n_cols
        for idx in M.colptr[j]:(M.colptr[j + 1] - 1)
            if !iszero(M.nzval[idx])
                algebraic_vars[j] = false
                algebraic_eqs[M.rowval[idx]] = false
            end
        end
    end

    return algebraic_vars, algebraic_eqs
end

# Fallback for non-sparse matrices (original behavior)
function find_algebraic_vars_eqs(M::AbstractMatrix)
    algebraic_vars = vec(all(iszero, M, dims = 1))
    algebraic_eqs = vec(all(iszero, M, dims = 2))
    return algebraic_vars, algebraic_eqs
end

# Handle SciMLOperators (e.g., MatrixOperator) by converting to matrix
function find_algebraic_vars_eqs(M::AbstractSciMLOperator)
    return find_algebraic_vars_eqs(convert(AbstractMatrix, M))
end

# Optimized tolerance checking that avoids allocations
@inline function check_dae_tolerance(integrator, err, abstol, t, ::Val{true})
    if abstol isa Number
        return integrator.opts.internalnorm(err, t) / abstol <= 1
    else
        @. err = err / abstol  # Safe for in-place functions
        return integrator.opts.internalnorm(err, t) <= 1
    end
end

@inline function check_dae_tolerance(integrator, err, abstol, t, ::Val{false})
    if abstol isa Number
        return integrator.opts.internalnorm(err, t) / abstol <= 1
    else
        return integrator.opts.internalnorm(err ./ abstol, t) <= 1  # Allocates for out-of-place
    end
end

function default_nlsolve(
        ::Nothing, isinplace::Val{true}, u, ::AbstractNonlinearProblem, autodiff = false)
    FastShortcutNonlinearPolyalg(;
        autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(
        ::Nothing, isinplace::Val{true}, u, ::NonlinearLeastSquaresProblem, autodiff = false)
    FastShortcutNLLSPolyalg(; autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(
        ::Nothing, isinplace::Val{false}, u, ::AbstractNonlinearProblem, autodiff = false)
    FastShortcutNonlinearPolyalg(;
        autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(
        ::Nothing, isinplace::Val{false}, u, ::NonlinearLeastSquaresProblem, autodiff = false)
    FastShortcutNLLSPolyalg(; autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(::Nothing, isinplace::Val{false}, u::StaticArray,
        ::AbstractNonlinearProblem, autodiff = false)
    SimpleTrustRegion(autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(::Nothing, isinplace::Val{false}, u::StaticArray,
        ::NonlinearLeastSquaresProblem, autodiff = false)
    SimpleGaussNewton(autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end

## ShampineCollocationInit

#=
The method:

du = (u-u0)/h
Solve for `u`

=#

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator,
        prob::ODEProblem, alg::DiffEqBase.ShampineCollocationInit,
        isinplace::Val{true})
    (; p, t, f) = integrator
    M = integrator.f.mass_matrix
    dtmax = integrator.opts.dtmax
    tmp = first(get_tmp_cache(integrator))
    u0 = integrator.u

    initdt = alg.initdt
    dt = if initdt === nothing
        integrator.dt != 0 ? min(integrator.dt / 5, dtmax) :
        (prob.tspan[end] - prob.tspan[begin]) / 1000 # Haven't implemented norm reduction
    else
        initdt
    end

    algebraic_vars, algebraic_eqs = find_algebraic_vars_eqs(M)
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    update_coefficients!(M, u0, p, t)
    f(tmp, u0, p, t)
    tmp .= ArrayInterface.restructure(tmp, algebraic_eqs .* _vec(tmp))

    check_dae_tolerance(integrator, tmp, integrator.opts.abstol, t, isinplace) && return

    if isdefined(integrator.cache, :nlsolver) && !isnothing(alg.nlsolve)
        # backward Euler
        nlsolver = integrator.cache.nlsolver
        oldγ, oldc, oldmethod,
        olddt = nlsolver.γ, nlsolver.c, nlsolver.method,
        integrator.dt
        nlsolver.tmp .= integrator.uprev
        nlsolver.γ, nlsolver.c = 1, 1
        nlsolver.method = DIRK
        integrator.dt = dt
        z = nlsolve!(nlsolver, integrator, integrator.cache)
        nlsolver.γ, nlsolver.c, nlsolver.method,
        integrator.dt = oldγ, oldc, oldmethod,
        olddt
        failed = nlsolvefail(nlsolver)
        @.. broadcast=false integrator.u=integrator.uprev+z
    else

        # _u0 should be non-dual since NonlinearSolve does not differentiate the solver
        # These non-dual values are thus used to make the caches
        #_du = SciMLBase.value.(du)
        _u0 = DiffEqBase.value.(u0)

        # If not doing auto-diff of the solver, save an allocation
        if typeof(u0) === typeof(_u0)
            tmp = get_tmp_cache(integrator)[1]
        else
            tmp = copy(_u0)
        end

        isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff ||
               typeof(u0) !== typeof(_u0)
        if isAD
            chunk = ForwardDiff.pickchunksize(length(tmp))
            _tmp = dualcache(tmp, chunk)
        else
            _tmp = tmp
        end

        nlequation! = @closure (out, u, p) -> begin
            TP = DiffEqBase.anyeltypedual(p)
            if TP <: Dual
                T = Base.promote_type(eltype(u), TP)
            else
                T = eltype(u)
            end
            update_coefficients!(M, u, p, t)
            # f(u,p,t) + M * (u0 - u)/dt
            tmp = isAD ? get_tmp(_tmp, T) : _tmp
            @. tmp = (_u0 - u) / dt
            mul!(_vec(out), M, _vec(tmp))
            f(tmp, u, p, t)
            out .+= tmp
            nothing
        end

        jac = if isnothing(f.jac)
            f.jac
        else
            @closure (J, u, p) -> begin
                # f(u,p,t) + M * (u0 - u)/dt
                # df(u,p,t)/du - M/dt
                f.jac(J, u, p, t)
                J .-= M .* inv(dt)
                nothing
            end
        end

        nlfunc = NonlinearFunction(nlequation!;
            jac_prototype = f.jac_prototype,
            jac = jac)
        nlprob = NonlinearProblem(nlfunc, integrator.u, p)
        nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0, nlprob, isAD)
        nlsol = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
            reltol = integrator.opts.reltol)
        integrator.u .= nlsol.u
        failed = nlsol.retcode != ReturnCode.Success
    end
    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
    end

    if failed
        @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator,
        prob::ODEProblem, alg::DiffEqBase.ShampineCollocationInit,
        isinplace::Val{false})
    (; p, t, f) = integrator
    u0 = integrator.u
    M = integrator.f.mass_matrix
    dtmax = integrator.opts.dtmax

    initdt = alg.initdt
    dt = if initdt === nothing
        integrator.dt != 0 ? min(integrator.dt / 5, dtmax) :
        (prob.tspan[end] - prob.tspan[begin]) / 1000 # Haven't implemented norm reduction
    else
        initdt
    end

    algebraic_vars, algebraic_eqs = find_algebraic_vars_eqs(M)
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    update_coefficients!(M, u0, p, t)
    du = f(u0, p, t)
    resid = _vec(du)[algebraic_eqs]

    check_dae_tolerance(integrator, resid, integrator.opts.abstol, t, isinplace) && return

    if isdefined(integrator.cache, :nlsolver) && !isnothing(alg.nlsolve)
        # backward Euler
        nlsolver = integrator.cache.nlsolver
        oldγ, oldc, oldmethod,
        olddt = nlsolver.γ, nlsolver.c, nlsolver.method,
        integrator.dt
        nlsolver.tmp .= integrator.uprev
        nlsolver.γ, nlsolver.c = 1, 1
        nlsolver.method = DIRK
        integrator.dt = dt
        z = nlsolve!(nlsolver, integrator, integrator.cache)
        nlsolver.γ, nlsolver.c, nlsolver.method,
        integrator.dt = oldγ, oldc, oldmethod,
        olddt
        failed = nlsolvefail(nlsolver)
        @.. broadcast=false integrator.u=integrator.uprev+z
    else
        nlequation_oop = @closure (u, _) -> begin
            update_coefficients!(M, u, p, t)
            M * (u - u0) / dt - f(u, p, t)
        end

        jac = if isnothing(f.jac)
            f.jac
        else
            @closure (u, p) -> begin
                return M * (u .- u0) ./ dt .- f.jac(u, p, t)
            end
        end

        nlfunc = NonlinearFunction(nlequation_oop;
            jac_prototype = f.jac_prototype,
            jac = jac)
        nlprob = NonlinearProblem(nlfunc, u0)
        nlsolve = default_nlsolve(alg.nlsolve, isinplace, nlprob, u0)

        nlsol = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
            reltol = integrator.opts.reltol)
        integrator.u = nlsol.u
        failed = nlsol.retcode != ReturnCode.Success
    end

    integrator.uprev = copy(integrator.u)
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = copy(integrator.uprev)
    end

    if failed
        @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator, prob::DAEProblem,
        alg::ShampineCollocationInit, isinplace::Val{true})
    (; p, t, f) = integrator
    u0 = integrator.u

    dtmax = integrator.opts.dtmax
    resid = get_tmp_cache(integrator)[2]

    dt = t != 0 ? min(t / 1000, dtmax / 10) : dtmax / 10 # Haven't implemented norm reduction

    f(resid, integrator.du, u0, p, t)
    check_dae_tolerance(integrator, resid, integrator.opts.abstol, t, isinplace) && return

    # _du and _u should be non-dual since NonlinearSolve does not differentiate the solver
    # These non-dual values are thus used to make the caches
    #_du = SciMLBase.value.(du)
    _u0 = DiffEqBase.value.(u0)

    # If not doing auto-diff of the solver, save an allocation
    if typeof(u0) === typeof(_u0)
        tmp = get_tmp_cache(integrator)[1]
    else
        tmp = copy(_u0)
    end

    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff || typeof(u0) !== typeof(_u0)
    if isAD
        chunk = ForwardDiff.pickchunksize(length(tmp))
        _tmp = dualcache(tmp, chunk)
    else
        _tmp = tmp
    end

    nlequation! = @closure (out, u, p) -> begin
        TP = DiffEqBase.anyeltypedual(p)
        if TP <: Dual
            T = Base.promote_type(eltype(u), TP)
        else
            T = eltype(u)
        end
        tmp = isAD ? get_tmp(_tmp, T) : _tmp
        #M * (u-u0)/dt - f(u,p,t)
        @. tmp = (u - _u0) / dt
        f(out, tmp, u, p, t)
        nothing
    end

    jac = if isnothing(f.jac)
        f.jac
    else
        @closure (J, u, p) -> begin
            f.jac(J, u, p, inv(dt), t)
            nothing
        end
    end

    nlfunc = NonlinearFunction(nlequation!;
        jac_prototype = f.jac_prototype,
        jac = jac)
    nlprob = NonlinearProblem(nlfunc, u0, p)
    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0, nlprob, isAD)
    nlsol = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
        reltol = integrator.opts.reltol)

    integrator.u = nlsol.u
    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
    end
    if nlsol.retcode != ReturnCode.Success
        @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator, prob::DAEProblem,
        alg::ShampineCollocationInit, isinplace::Val{false})
    (; p, t, f) = integrator
    u0 = integrator.u
    dtmax = integrator.opts.dtmax

    dt = t != 0 ? min(t / 1000, dtmax / 10) : dtmax / 10 # Haven't implemented norm reduction

    nlequation_oop = u -> begin
        f((u - u0) / dt, u, p, t)
    end

    nlequation = (u, _) -> nlequation_oop(u)

    resid = f(integrator.du, u0, p, t)
    check_dae_tolerance(integrator, resid, integrator.opts.abstol, t, isinplace) && return

    jac = if isnothing(f.jac)
        f.jac
    else
        @closure (u, p) -> begin
            return f.jac(u, p, inv(dt), t)
        end
    end
    nlfunc = NonlinearFunction(nlequation; jac_prototype = f.jac_prototype,
        jac = jac)
    nlprob = NonlinearProblem(nlfunc, u0)
    nlsolve = default_nlsolve(alg.nlsolve, isinplace, nlprob, u0)

    nlfunc = NonlinearFunction(nlequation; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, u0)
    nlsol = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
        reltol = integrator.opts.reltol)

    integrator.u = nlsol.u

    integrator.uprev = copy(integrator.u)
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = copy(integrator.uprev)
    end
    if nlsol.retcode != ReturnCode.Success
        @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

## BrownFullBasic

#=
The method:

Keep differential variables constant
Solve for the algebraic variables

=#

algebraic_jacobian(::Nothing, algebraic_eqs, algebraic_vars) = nothing
function algebraic_jacobian(jac_prototype::T, algebraic_eqs,
        algebraic_vars) where {T <: AbstractMatrix}
    jac_prototype[algebraic_eqs, algebraic_vars]
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator, prob::ODEProblem,
        alg::DiffEqBase.BrownFullBasicInit, isinplace::Val{true})
    (; p, t, f) = integrator
    u = integrator.u
    M = integrator.f.mass_matrix
    M isa UniformScaling && return
    update_coefficients!(M, u, p, t)
    algebraic_vars, algebraic_eqs = find_algebraic_vars_eqs(M)

    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    tmp = get_tmp_cache(integrator)[1]

    f(tmp, u, p, t)

    tmp .= ArrayInterface.restructure(tmp, algebraic_eqs .* _vec(tmp))

    check_dae_tolerance(integrator, tmp, alg.abstol, t, isinplace) && return
    alg_u = @view u[algebraic_vars]

    # These non-dual values are thus used to make the caches
    _u = DiffEqBase.value.(u)

    # If auto-diff of the solver, should be non-dual since NonlinearSolve does not differentiate the solver
    if typeof(u) !== typeof(_u)
        tmp = DiffEqBase.value.(tmp)
    end

    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff || typeof(u) !== typeof(_u)
    if isAD
        csize = count(algebraic_vars)
        if !(p isa SciMLBase.NullParameters) && typeof(_u) !== typeof(u)
            if isscimlstructure(p)
                csize = max(csize, length(canonicalize(Tunable(), p)[1]))
            else
                csize = max(csize, length(p))
            end
        end
        chunk = ForwardDiff.pickchunksize(csize)
        _tmp = dualcache(tmp, chunk)
        _du_tmp = dualcache(similar(tmp), chunk)
    else
        _tmp, _du_tmp = tmp, similar(tmp)
    end

    nlequation! = @closure (out, x, p) -> begin
        TP = DiffEqBase.anyeltypedual(p)
        if TP <: Dual
            T = Base.promote_type(eltype(x), TP)
        else
            T = eltype(x)
        end
        uu = isAD ? get_tmp(_tmp, T) : _tmp
        du_tmp = isAD ? get_tmp(_du_tmp, T) : _du_tmp
        copyto!(uu, _u)
        alg_uu = @view uu[algebraic_vars]
        alg_uu .= x
        f(du_tmp, uu, p, t)
        out .= @view du_tmp[algebraic_eqs]
        return nothing
    end

    J = algebraic_jacobian(f.jac_prototype, algebraic_eqs, algebraic_vars)
    nlfunc = NonlinearFunction(nlequation!; jac_prototype = J)
    nlprob = NonlinearProblem(nlfunc, alg_u, p)
    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u, nlprob, isAD)

    nlsol = solve(nlprob, nlsolve; abstol = alg.abstol, reltol = integrator.opts.reltol)
    alg_u .= nlsol.u

    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
    end

    if nlsol.retcode != ReturnCode.Success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator, prob::ODEProblem,
        alg::DiffEqBase.BrownFullBasicInit, isinplace::Val{false})
    (; p, t, f) = integrator

    u0 = integrator.u
    M = integrator.f.mass_matrix
    update_coefficients!(M, u0, p, t)
    algebraic_vars, algebraic_eqs = find_algebraic_vars_eqs(M)
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return

    du = f(u0, p, t)
    resid = _vec(du)[algebraic_eqs]

    check_dae_tolerance(integrator, resid, alg.abstol, t, isinplace) && return

    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff
    if isAD
        chunk = ForwardDiff.pickchunksize(count(algebraic_vars))
        _tmp = dualcache(similar(u0), chunk)
    else
        _tmp = similar(u0)
    end

    if u0 isa Number
        # This doesn't fix static arrays!
        u = [u0]
    else
        u = u0
    end

    nlequation = @closure (x, _) -> begin
        uu = isAD ? get_tmp(_tmp, x) : _tmp
        copyto!(uu, integrator.u)
        alg_u = @view uu[algebraic_vars]
        alg_u .= x
        du = f(uu, p, t)
        @views du[algebraic_eqs]
    end

    J = algebraic_jacobian(f.jac_prototype, algebraic_eqs, algebraic_vars)
    nlfunc = NonlinearFunction(nlequation; jac_prototype = J)
    nlprob = NonlinearProblem(nlfunc, u0[algebraic_vars])
    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0, nlprob, isAD)

    nlsol = solve(nlprob, nlsolve)

    u[algebraic_vars] .= nlsol.u

    if u0 isa Number
        # This doesn't fix static arrays!
        integrator.u = first(u)
    else
        integrator.u = u
    end

    integrator.uprev = copy(integrator.u)
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = copy(integrator.uprev)
    end

    if nlsol.retcode != ReturnCode.Success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator, prob::DAEProblem,
        alg::DiffEqBase.BrownFullBasicInit, isinplace::Val{true})
    (; p, t, f) = integrator
    differential_vars = prob.differential_vars
    u = integrator.u
    du = integrator.du

    # _du and _u should be non-dual since NonlinearSolve does not differentiate the solver
    # These non-dual values are thus used to make the caches
    _du = DiffEqBase.value.(du)
    _u = DiffEqBase.value.(u)

    # If not doing auto-diff of the solver, save an allocation
    if typeof(u) === typeof(_u)
        tmp = get_tmp_cache(integrator)[1]
        du_tmp = get_tmp_cache(integrator)[2]
    else
        tmp = copy(_u)
        du_tmp = copy(_du)
    end

    # Can be the same as tmp
    normtmp = get_tmp_cache(integrator)[1]
    f(normtmp, du, u, p, t)

    if check_dae_tolerance(integrator, normtmp, alg.abstol, t, isinplace)
        return
    elseif differential_vars === nothing
        error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions, differential_vars, or use a different initialization algorithm.")
    end

    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff || typeof(u) !== typeof(_u)
    if isAD
        chunk = ForwardDiff.pickchunksize(length(tmp))
        _tmp = dualcache(tmp, chunk)
        _du_tmp = dualcache(du_tmp, chunk)
    else
        _tmp, _du_tmp = tmp, du_tmp
    end

    nlequation! = @closure (out, x, p) -> begin
        TP = DiffEqBase.anyeltypedual(p)
        if TP <: Dual
            T = Base.promote_type(eltype(x), TP)
        else
            T = eltype(x)
        end
        du_tmp = isAD ? get_tmp(_du_tmp, T) : _du_tmp
        uu = isAD ? get_tmp(_tmp, T) : _tmp

        @. du_tmp = ifelse(differential_vars, x, _du)
        @. uu = ifelse(differential_vars, _u, x)

        f(out, du_tmp, uu, p, t)
    end

    nlsolve_alg = alg.nlsolve
    if nlsolve_alg !== nothing
        nlsolve = nlsolve_alg
    else
        nlsolve = NewtonRaphson(autodiff = alg_autodiff(integrator.alg))
    end

    nlfunc = NonlinearFunction(nlequation!; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, ifelse.(differential_vars, du, u), p)
    nlsol = solve(nlprob, nlsolve; abstol = alg.abstol, reltol = integrator.opts.reltol)

    @. du = ifelse(differential_vars, nlsol.u, du)
    @. u = ifelse(differential_vars, u, nlsol.u)

    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
    end

    if nlsol.retcode != ReturnCode.Success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end

function _initialize_dae!(integrator::OrdinaryDiffEqCore.ODEIntegrator, prob::DAEProblem,
        alg::DiffEqBase.BrownFullBasicInit, isinplace::Val{false})
    (; p, t, f) = integrator
    differential_vars = prob.differential_vars

    if check_dae_tolerance(
        integrator, f(integrator.du, integrator.u, p, t), alg.abstol, t, isinplace)
        return
    elseif differential_vars === nothing
        error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions, differential_vars, or use a different initialization algorithm.")
    end

    if integrator.u isa Number && integrator.du isa Number
        # This doesn't fix static arrays!
        u = [integrator.u]
        du = [integrator.du]
    else
        u = integrator.u
        du = integrator.du
    end

    nlequation = @closure (x, _) -> begin
        du_ = ifelse.(differential_vars, x, du)
        u_ = ifelse.(differential_vars, u, x)
        f.f(du_, u_, p, t)
    end

    nlfunc = NonlinearFunction(nlequation; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, ifelse.(differential_vars, du, u))

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, nlprob, integrator.u)

    @show nlsolve

    nlsol = solve(nlprob, nlsolve)

    du = ifelse.(differential_vars, nlsol.u, du)
    u = ifelse.(differential_vars, u, nlsol.u)

    if integrator.u isa Number && integrator.du isa Number
        # This doesn't fix static arrays!
        integrator.u = first(u)
        integrator.du = first(du)
    else
        integrator.u = u
        integrator.du = du
    end

    integrator.uprev = copy(integrator.u)
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = copy(integrator.uprev)
    end

    if nlsol.retcode != ReturnCode.Success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
    return
end
