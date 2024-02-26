struct DefaultInit <: DiffEqBase.DAEInitializationAlgorithm end

struct ShampineCollocationInit{T, F} <: DiffEqBase.DAEInitializationAlgorithm
    initdt::T
    nlsolve::F
end
function ShampineCollocationInit(; initdt = nothing, nlsolve = nothing)
    ShampineCollocationInit(initdt, nlsolve)
end
function ShampineCollocationInit(initdt)
    ShampineCollocationInit(; initdt = initdt, nlsolve = nothing)
end

struct BrownFullBasicInit{T, F} <: DiffEqBase.DAEInitializationAlgorithm
    abstol::T
    nlsolve::F
end
function BrownFullBasicInit(; abstol = 1e-10, nlsolve = nothing)
    BrownFullBasicInit(abstol, nlsolve)
end
BrownFullBasicInit(abstol) = BrownFullBasicInit(; abstol = abstol, nlsolve = nothing)

default_nlsolve(alg, isinplace, u, autodiff = false) = alg
function default_nlsolve(::Nothing, isinplace, u, ::NonlinearProblem, autodiff = false)
    FastShortcutNonlinearPolyalg(;
        autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(::Nothing, isinplace::Val{false}, u::StaticArray,
        ::NonlinearProblem, autodiff = false)
    SimpleTrustRegion(autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end

function default_nlsolve(
        ::Nothing, isinplace, u, ::NonlinearLeastSquaresProblem, autodiff = false)
    FastShortcutNLLSPolyalg(; autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end
function default_nlsolve(::Nothing, isinplace::Val{false}, u::StaticArray,
        ::NonlinearLeastSquaresProblem, autodiff = false)
    SimpleGaussNewton(autodiff = autodiff ? AutoForwardDiff() : AutoFiniteDiff())
end

struct OverrideInit{T, F} <: DiffEqBase.DAEInitializationAlgorithm
    abstol::T
    nlsolve::F
end

function OverrideInit(; abstol = 1e-10, nlsolve = nothing)
    OverrideInit(abstol, nlsolve)
end
OverrideInit(abstol) = OverrideInit(; abstol = abstol, nlsolve = nothing)

## Notes

#=
differential_vars = [any(!iszero,x) for x in eachcol(M)]

A column should be zero for an algebraic variable, since that means that the
derivative term doesn't show up in any equations (i.e. is an algebraic variable).
The rows are not necessarily non-zero, for example a flux condition between two
differential variables. But if it's a condition that doesn't involve the algebraic
variable, then the system is not Index 1!

=#

## Expansion

function DiffEqBase.initialize_dae!(integrator::ODEIntegrator,
        initializealg = integrator.initializealg)
    _initialize_dae!(integrator, integrator.sol.prob,
        initializealg,
        Val(DiffEqBase.isinplace(integrator.sol.prob)))
end

## Default algorithms

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::DefaultInit, x::Val{true})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::DefaultInit, x::Val{false})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::DefaultInit, x::Val{false})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    elseif prob.differential_vars === nothing
        _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::DefaultInit, x::Val{true})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    elseif prob.differential_vars === nothing
        _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

## NoInit

function _initialize_dae!(integrator, prob::Union{ODEProblem, DAEProblem},
        alg::NoInit, x::Union{Val{true}, Val{false}})
end

## OverrideInit

function _initialize_dae!(integrator, prob::Union{ODEProblem, DAEProblem},
        alg::OverrideInit, isinplace::Union{Val{true}, Val{false}})
    initializeprob = prob.f.initializeprob
    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff
    alg = default_nlsolve(alg.nlsolve, isinplace, initializeprob.u0, initializeprob, isAD)
    nlsol = solve(initializeprob, alg)
    if isinplace === Val{true}()
        integrator.u .= prob.f.initializeprobmap(nlsol)
    elseif isinplace === Val{false}()
        integrator.u = prob.f.initializeprobmap(nlsol)
    else
        error("Unreachable reached. Report this error.")
    end

    if nlsol.retcode != ReturnCode.Success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
end

## ShampineCollocationInit

#=
The method:

du = (u-u0)/h
Solve for `u`

=#

function _initialize_dae!(integrator, prob::ODEProblem, alg::ShampineCollocationInit,
        isinplace::Val{true})
    @unpack p, t, f = integrator
    M = integrator.f.mass_matrix
    dtmax = integrator.opts.dtmax
    tmp = first(get_tmp_cache(integrator))
    u0 = integrator.u

    dt = if alg.initdt === nothing
        integrator.dt != 0 ? min(integrator.dt / 5, dtmax) :
        (prob.tspan[end] - prob.tspan[begin]) / 1000 # Haven't implemented norm reduction
    else
        alg.initdt
    end

    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    update_coefficients!(M, u0, p, t)
    f(tmp, u0, p, t)
    tmp .= ArrayInterface.restructure(tmp, algebraic_eqs .* _vec(tmp))

    integrator.opts.internalnorm(tmp, t) <= integrator.opts.abstol && return

    if isdefined(integrator.cache, :nlsolver) && !isnothing(alg.nlsolve)
        # backward Euler
        nlsolver = integrator.cache.nlsolver
        oldγ, oldc, oldmethod, olddt = nlsolver.γ, nlsolver.c, nlsolver.method,
        integrator.dt
        nlsolver.tmp .= integrator.uprev
        nlsolver.γ, nlsolver.c = 1, 1
        nlsolver.method = DIRK
        integrator.dt = dt
        z = nlsolve!(nlsolver, integrator, integrator.cache)
        nlsolver.γ, nlsolver.c, nlsolver.method, integrator.dt = oldγ, oldc, oldmethod,
        olddt
        failed = nlsolvefail(nlsolver)
        @.. broadcast=false integrator.u=integrator.uprev + z
    else

        # _u0 should be non-dual since NonlinearSolve does not differentiate the solver
        # These non-dual values are thus used to make the caches
        #_du = DiffEqBase.value.(du)
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
            _tmp = PreallocationTools.dualcache(tmp, chunk)
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
            tmp = isAD ? PreallocationTools.get_tmp(_tmp, T) : _tmp
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

function _initialize_dae!(integrator, prob::ODEProblem, alg::ShampineCollocationInit,
        isinplace::Val{false})
    @unpack p, t, f = integrator
    u0 = integrator.u
    M = integrator.f.mass_matrix
    dtmax = integrator.opts.dtmax

    dt = if alg.initdt === nothing
        integrator.dt != 0 ? min(integrator.dt / 5, dtmax) :
        (prob.tspan[end] - prob.tspan[begin]) / 1000 # Haven't implemented norm reduction
    else
        alg.initdt
    end

    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    update_coefficients!(M, u0, p, t)
    du = f(u0, p, t)
    resid = _vec(du)[algebraic_eqs]

    integrator.opts.internalnorm(resid, t) <= integrator.opts.abstol && return

    if isdefined(integrator.cache, :nlsolver) && !isnothing(alg.nlsolve)
        # backward Euler
        nlsolver = integrator.cache.nlsolver
        oldγ, oldc, oldmethod, olddt = nlsolver.γ, nlsolver.c, nlsolver.method,
        integrator.dt
        nlsolver.tmp .= integrator.uprev
        nlsolver.γ, nlsolver.c = 1, 1
        nlsolver.method = DIRK
        integrator.dt = dt
        z = nlsolve!(nlsolver, integrator, integrator.cache)
        nlsolver.γ, nlsolver.c, nlsolver.method, integrator.dt = oldγ, oldc, oldmethod,
        olddt
        failed = nlsolvefail(nlsolver)
        @.. broadcast=false integrator.u=integrator.uprev + z
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

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::ShampineCollocationInit, isinplace::Val{true})
    @unpack p, t, f = integrator
    u0 = integrator.u

    dtmax = integrator.opts.dtmax
    resid = get_tmp_cache(integrator)[2]

    dt = t != 0 ? min(t / 1000, dtmax / 10) : dtmax / 10 # Haven't implemented norm reduction

    f(resid, integrator.du, u0, p, t)
    integrator.opts.internalnorm(resid, t) <= integrator.opts.abstol && return

    # _du and _u should be non-dual since NonlinearSolve does not differentiate the solver
    # These non-dual values are thus used to make the caches
    #_du = DiffEqBase.value.(du)
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
        _tmp = PreallocationTools.dualcache(tmp, chunk)
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
        tmp = isAD ? PreallocationTools.get_tmp(_tmp, T) : _tmp
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

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::ShampineCollocationInit, isinplace::Val{false})
    @unpack p, t, f = integrator
    u0 = integrator.u
    dtmax = integrator.opts.dtmax

    dt = t != 0 ? min(t / 1000, dtmax / 10) : dtmax / 10 # Haven't implemented norm reduction

    nlequation_oop = u -> begin
        f((u - u0) / dt, u, p, t)
    end

    nlequation = (u, _) -> nlequation_oop(u)

    resid = f(integrator.du, u0, p, t)
    integrator.opts.internalnorm(resid, t) <= integrator.opts.abstol && return

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

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::BrownFullBasicInit, isinplace::Val{true})
    @unpack p, t, f = integrator
    u = integrator.u
    M = integrator.f.mass_matrix
    M isa UniformScaling && return
    update_coefficients!(M, u, p, t)
    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    tmp = get_tmp_cache(integrator)[1]

    f(tmp, u, p, t)

    tmp .= ArrayInterface.restructure(tmp, algebraic_eqs .* _vec(tmp))

    integrator.opts.internalnorm(tmp, t) <= alg.abstol && return
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
            try
                csize = max(csize, length(p))
            catch
            end
        end
        chunk = ForwardDiff.pickchunksize(csize)
        _tmp = PreallocationTools.dualcache(tmp, chunk)
        _du_tmp = PreallocationTools.dualcache(similar(tmp), chunk)
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
        uu = isAD ? PreallocationTools.get_tmp(_tmp, T) : _tmp
        du_tmp = isAD ? PreallocationTools.get_tmp(_du_tmp, T) : _du_tmp
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
    alg_u .= nlsol

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

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::BrownFullBasicInit, isinplace::Val{false})
    @unpack p, t, f = integrator

    u0 = integrator.u
    M = integrator.f.mass_matrix
    update_coefficients!(M, u0, p, t)
    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return

    du = f(u0, p, t)
    resid = _vec(du)[algebraic_eqs]

    integrator.opts.internalnorm(resid, t) <= alg.abstol && return

    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff
    if isAD
        chunk = ForwardDiff.pickchunksize(count(algebraic_vars))
        _tmp = PreallocationTools.dualcache(similar(u0), chunk)
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
        uu = isAD ? PreallocationTools.get_tmp(_tmp, x) : _tmp
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

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::BrownFullBasicInit, isinplace::Val{true})
    @unpack p, t, f = integrator
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

    if integrator.opts.internalnorm(normtmp, t) <= alg.abstol
        return
    elseif differential_vars === nothing
        error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions, differential_vars, or use a different initialization algorithm.")
    end

    isAD = alg_autodiff(integrator.alg) isa AutoForwardDiff || typeof(u) !== typeof(_u)
    if isAD
        chunk = ForwardDiff.pickchunksize(length(tmp))
        _tmp = PreallocationTools.dualcache(tmp, chunk)
        _du_tmp = PreallocationTools.dualcache(du_tmp, chunk)
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
        du_tmp = isAD ? PreallocationTools.get_tmp(_du_tmp, T) : _du_tmp
        uu = isAD ? PreallocationTools.get_tmp(_tmp, T) : _tmp

        @. du_tmp = ifelse(differential_vars, x, _du)
        @. uu = ifelse(differential_vars, _u, x)

        f(out, du_tmp, uu, p, t)
    end

    if alg.nlsolve !== nothing
        nlsolve = alg.nlsolve
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

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::BrownFullBasicInit, isinplace::Val{false})
    @unpack p, t, f = integrator
    differential_vars = prob.differential_vars

    if integrator.opts.internalnorm(f(integrator.du, integrator.u, p, t), t) <= alg.abstol
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
        du = ifelse.(differential_vars, x, du)
        u = ifelse.(differential_vars, u, x)
        f(du, u, p, t)
    end

    nlfunc = NonlinearFunction(nlequation; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, ifelse.(differential_vars, du, u))

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, nlprob, integrator.u)

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
