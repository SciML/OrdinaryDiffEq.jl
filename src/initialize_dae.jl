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

using SciMLNLSolve
default_nlsolve(alg, isinplace, u, autodiff = false) = alg
function default_nlsolve(::Nothing, isinplace, u, autodiff = false)
    NLSolveJL(autodiff = autodiff)
end

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
    _initialize_dae!(integrator, prob,
                     BrownFullBasicInit(integrator.opts.abstol), x)
end

function _initialize_dae!(integrator, prob::ODEProblem,
                          alg::DefaultInit, x::Val{false})
    _initialize_dae!(integrator, prob,
                     BrownFullBasicInit(integrator.opts.abstol), x)
end

function _initialize_dae!(integrator, prob::DAEProblem,
                          alg::DefaultInit, x::Val{false})
    if prob.differential_vars === nothing
        _initialize_dae!(integrator, prob,
                         ShampineCollocationInit(), x)
    else
        _initialize_dae!(integrator, prob,
                         BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
                          alg::DefaultInit, x::Val{true})
    if prob.differential_vars === nothing
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
        integrator.dt != 0 ? min(integrator.dt / 5, dtmax) : 1 // 1000 # Haven't implemented norm reduction
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

    if isdefined(integrator.cache, :nlsolver)
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
        # TODO: failure handling
        nlsolvefail(nlsolver) &&
            @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
        @.. broadcast=false integrator.u=integrator.uprev + z
    else
        isad = alg_autodiff(integrator.alg)
        if isad
            chunk = ForwardDiff.pickchunksize(length(tmp))
            _tmp = PreallocationTools.dualcache(tmp, chunk)
        else
            _tmp = tmp
        end

        nlequation! = @closure (out, u, p) -> begin
            update_coefficients!(M, u, p, t)
            #M * (u-u0)/dt - f(u,p,t)
            tmp = isad ? PreallocationTools.get_tmp(_tmp, u) : _tmp
            @. tmp = (u - u0) / dt
            mul!(_vec(out), M, _vec(tmp))
            f(tmp, u, p, t)
            out .-= tmp
            nothing
        end

        nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0, isad)

        nlfunc = NonlinearFunction(nlequation!;
                                   jac_prototype = f.jac_prototype)
        nlprob = NonlinearProblem(nlfunc, integrator.u, p)
        r = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
                  reltol = integrator.opts.reltol)
        integrator.u .= r.u
    end
    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
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
        integrator.dt != 0 ? min(integrator.dt / 5, dtmax) : 1 // 1000 # Haven't implemented norm reduction
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

    if isdefined(integrator.cache, :nlsolver)
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
        # TODO: failure handling
        nlsolvefail(nlsolver) &&
            @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
        @.. broadcast=false integrator.u=integrator.uprev + z
    else
        nlequation_oop = @closure (u, _) -> begin
            update_coefficients!(M, u, p, t)
            M * (u - u0) / dt - f(u, p, t)
        end

        nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0)

        nlfunc = NonlinearFunction(nlequation_oop;
                                   jac_prototype = f.jac_prototype)
        nlprob = NonlinearProblem(nlfunc, u0)
        nlsol = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
                      reltol = integrator.opts.reltol)
        integrator.u = nlsol.u
    end
    integrator.uprev = integrator.u
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = integrator.uprev
    end

    return
end

function _initialize_dae!(integrator, prob::DAEProblem,
                          alg::ShampineCollocationInit, isinplace::Val{true})
    @unpack p, t, f = integrator
    u0 = integrator.u

    dtmax = integrator.opts.dtmax
    tmp = get_tmp_cache(integrator)[1]
    resid = get_tmp_cache(integrator)[2]

    dt = t != 0 ? min(t / 1000, dtmax) : dtmax # Haven't implemented norm reduction

    f(resid, integrator.du, u0, p, t)
    integrator.opts.internalnorm(resid, t) <= integrator.opts.abstol && return

    isad = alg_autodiff(integrator.alg)
    if isad
        chunk = ForwardDiff.pickchunksize(length(tmp))
        _tmp = PreallocationTools.dualcache(tmp, chunk)
    else
        _tmp = tmp
    end

    nlequation! = @closure (out, u, p) -> begin
        tmp = isad ? PreallocationTools.get_tmp(_tmp, u) : _tmp
        #M * (u-u0)/dt - f(u,p,t)
        @. tmp = (u - u0) / dt
        f(out, tmp, u, p, t)
        nothing
    end

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0, isad)

    nlfunc = NonlinearFunction(nlequation!; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, u0, p)
    r = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
              reltol = integrator.opts.reltol)

    integrator.u = r.u
    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
    end

    return
end

function _initialize_dae!(integrator, prob::DAEProblem,
                          alg::ShampineCollocationInit, isinplace::Val{false})
    @unpack p, t, f = integrator
    u0 = integrator.u
    dtmax = integrator.opts.dtmax

    dt = t != 0 ? min(t / 1000, dtmax / 10) : dtmax # Haven't implemented norm reduction

    nlequation_oop = u -> begin f((u - u0) / dt, u, p, t) end

    nlequation = (u, _) -> nlequation_oop(u)

    resid = f(integrator.du, u0, p, t)
    integrator.opts.internalnorm(resid, t) <= integrator.opts.abstol && return

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0)

    nlfunc = NonlinearFunction(nlequation; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, u0)
    r = solve(nlprob, nlsolve; abstol = integrator.opts.abstol,
              reltol = integrator.opts.reltol)

    integrator.u = r.u
    integrator.uprev = integrator.u
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = integrator.uprev
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
    update_coefficients!(M, u, p, t)
    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    tmp = get_tmp_cache(integrator)[1]

    f(tmp, u, p, t)

    tmp .= ArrayInterface.restructure(tmp, algebraic_eqs .* _vec(tmp))

    integrator.opts.internalnorm(tmp, t) <= alg.abstol && return
    alg_u = @view u[algebraic_vars]

    isad = alg_autodiff(integrator.alg)
    if isad
        chunk = ForwardDiff.pickchunksize(count(algebraic_vars))
        _tmp = PreallocationTools.dualcache(tmp, chunk)
        _du_tmp = PreallocationTools.dualcache(tmp, chunk)
    else
        _tmp, _du_tmp = tmp, similar(tmp)
    end

    nlequation! = @closure (out, x, p) -> begin
        uu = isad ? PreallocationTools.get_tmp(_tmp, x) : _tmp
        du_tmp = isad ? PreallocationTools.get_tmp(_du_tmp, x) : _du_tmp
        copyto!(uu, integrator.u)
        alg_uu = @view uu[algebraic_vars]
        alg_uu .= x
        f(du_tmp, uu, p, t)
        out .= @view du_tmp[algebraic_eqs]
        return nothing
    end

    J = algebraic_jacobian(f.jac_prototype, algebraic_eqs, algebraic_vars)

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u, isad)

    nlfunc = NonlinearFunction(nlequation!; jac_prototype = J)
    nlprob = NonlinearProblem(nlfunc, u[algebraic_vars], p)
    r = solve(nlprob, nlsolve; abstol = alg.abstol, reltol = integrator.opts.reltol)

    alg_u .= r

    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
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

    isad = alg_autodiff(integrator.alg)
    if isad
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
        uu = isad ? PreallocationTools.get_tmp(_tmp, x) : _tmp
        copyto!(uu, integrator.u)
        alg_u = @view uu[algebraic_vars]
        alg_u .= x
        du = f(uu, p, t)
        @views du[algebraic_eqs]
    end

    J = algebraic_jacobian(f.jac_prototype, algebraic_eqs, algebraic_vars)

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, u0, isad)

    nlfunc = NonlinearFunction(nlequation; jac_prototype = J)
    nlprob = NonlinearProblem(nlfunc, u0[algebraic_vars])
    r = solve(nlprob, nlsolve)

    u[algebraic_vars] .= r.u

    if u0 isa Number
        # This doesn't fix static arrays!
        integrator.u = first(u)
    else
        integrator.u = u
    end

    integrator.uprev = integrator.u
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = integrator.uprev
    end

    return
end

function _initialize_dae!(integrator, prob::DAEProblem,
                          alg::BrownFullBasicInit, isinplace::Val{true})
    @unpack p, t, f = integrator
    differential_vars = prob.differential_vars
    u = integrator.u
    du = integrator.du

    tmp = get_tmp_cache(integrator)[1]
    du_tmp = get_tmp_cache(integrator)[2]
    f(tmp, du, u, p, t)

    if integrator.opts.internalnorm(tmp, t) <= alg.abstol
        return
    elseif differential_vars === nothing
        error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions, differential_vars, or use a different initialization algorithm.")
    end

    isad = alg_autodiff(integrator.alg)
    if isad
        chunk = ForwardDiff.pickchunksize(length(tmp))
        _tmp = PreallocationTools.dualcache(tmp, chunk)
        _du_tmp = PreallocationTools.dualcache(du_tmp, chunk)
    else
        _tmp, _du_tmp = tmp, similar(tmp)
    end

    nlequation! = @closure (out, x, p) -> begin
        du_tmp = isad ? PreallocationTools.get_tmp(_du_tmp, x) : _du_tmp
        uu = isad ? PreallocationTools.get_tmp(_tmp, x) : _tmp

        @. du_tmp = ifelse(differential_vars, x, du)
        @. uu = ifelse(differential_vars, u, x)

        f(out, du_tmp, uu, p, t)
    end

    if alg.nlsolve !== nothing
        nlsolve = alg.nlsolve
    else
        nlsolve = NewtonRaphson(autodiff = isad)
    end

    nlfunc = NonlinearFunction(nlequation!; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, ifelse.(differential_vars, du, u), p)
    r = solve(nlprob, nlsolve; abstol = alg.abstol, reltol = integrator.opts.reltol)

    @. du = ifelse(differential_vars, r.u, du)
    @. u = ifelse(differential_vars, u, r.u)

    recursivecopy!(integrator.uprev, integrator.u)
    if alg_extrapolates(integrator.alg)
        recursivecopy!(integrator.uprev2, integrator.uprev)
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

    nlsolve = default_nlsolve(alg.nlsolve, isinplace, integrator.u)

    nlfunc = NonlinearFunction(nlequation; jac_prototype = f.jac_prototype)
    nlprob = NonlinearProblem(nlfunc, ifelse.(differential_vars, du, u))

    r = solve(nlprob, nlsolve)

    du = ifelse.(differential_vars, r.u, du)
    u = ifelse.(differential_vars, u, r.u)

    if integrator.u isa Number && integrator.du isa Number
        # This doesn't fix static arrays!
        integrator.u = first(u)
        integrator.du = first(du)
    else
        integrator.u = u
        integrator.du = du
    end

    integrator.uprev = integrator.u
    if alg_extrapolates(integrator.alg)
        integrator.uprev2 = integrator.uprev
    end

    return
end
