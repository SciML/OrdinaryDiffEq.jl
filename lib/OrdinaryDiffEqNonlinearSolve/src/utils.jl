# Get mass_matrix from function, defaulting to I if field doesn't exist
# This allows nlsolver to work with DiscreteFunction which lacks mass_matrix
get_mass_matrix(f) = hasproperty(f, :mass_matrix) ? f.mass_matrix : I

get_status(nlsolver::AbstractNLSolver) = nlsolver.status
get_new_W_γdt_cutoff(nlsolver::AbstractNLSolver) = nlsolver.cache.new_W_γdt_cutoff
# handle FIRK
get_new_W_γdt_cutoff(alg::NewtonAlgorithm) = alg.new_W_γdt_cutoff

"""
    nlsolvefail(nlsolver) -> Bool
    nlsolvefail(status::NLStatus) -> Bool

Return whether a nonlinear solve failed, i.e. whether the solver's
[`NLStatus`](@ref) is non-positive (`SlowConvergence` or worse).
"""
nlsolvefail(nlsolver::AbstractNLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(status::NLStatus) = Int8(status) <= 0

relax(nlsolver::AbstractNLSolver) = relax(nlsolver.alg)
relax(alg::NLNewton) = alg.relax
relax(_) = 0 // 1

isnewton(nlsolver::AbstractNLSolver) = isnewton(nlsolver.cache)
isnewton(::AbstractNLSolverCache) = false
isnewton(::Union{NLNewtonCache, NLNewtonConstantCache}) = true

is_always_new(nlsolver::AbstractNLSolver) = is_always_new(nlsolver.alg)
check_div(nlsolver::AbstractNLSolver) = check_div(nlsolver.alg)
check_div(alg) = isdefined(alg, :check_div) ? alg.check_div : true

isJcurrent(nlsolver::AbstractNLSolver, integrator) = integrator.t == nlsolver.cache.J_t
isfirstcall(nlsolver::AbstractNLSolver) = nlsolver.cache.firstcall
isfirststage(nlsolver::AbstractNLSolver) = nlsolver.cache.firststage
setfirststage!(nlsolver::AbstractNLSolver, val::Bool) = setfirststage!(nlsolver.cache, val)
function setfirststage!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}, val::Bool)
    return (nlcache.firststage = val)
end
setfirststage!(::Any, val::Bool) = nothing
"""
    markfirststage!(nlsolver)

Mark the nonlinear solver as being on the first implicit stage of the current
step (sets the cache's `firststage` flag to `true`). Used so predictor/`W`-reuse
logic can distinguish the first stage from later ones.
"""
markfirststage!(nlsolver::AbstractNLSolver) = setfirststage!(nlsolver, true)

getnfails(_) = 0
getnfails(nlsolver::AbstractNLSolver) = nlsolver.nfails

set_new_W!(nlsolver::AbstractNLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}, val::Bool)::Bool = nlcache.new_W = val
get_new_W!(nlsolver::AbstractNLSolver)::Bool = get_new_W!(nlsolver.cache)
get_new_W!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache})::Bool = nlcache.new_W
get_new_W!(::AbstractNLSolverCache)::Bool = true
get_new_W!(nlcache::NonlinearSolveCache)::Bool = nlcache.new_W

get_W(nlsolver::AbstractNLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}) = nlcache.W

set_W_γdt!(nlsolver::AbstractNLSolver, W_γdt) = set_W_γdt!(nlsolver.cache, W_γdt)
function set_W_γdt!(nlcache::Union{NLNewtonCache, NLNewtonConstantCache}, W_γdt)
    nlcache.W_γdt = W_γdt
    return W_γdt
end

du_cache(nlsolver::AbstractNLSolver) = du_cache(nlsolver.cache)
du_cache(::AbstractNLSolverCache) = nothing
du_cache(nlcache::Union{NLFunctionalCache, NLAndersonCache, NLNewtonCache}) = (nlcache.k,)

"""
    du_alias_or_new(nlsolver, rate_prototype)

Return a derivative buffer for the nonlinear solve: reuse the solver cache's
existing `du` buffer when it has one, otherwise allocate `zero(rate_prototype)`.
"""
function du_alias_or_new(nlsolver::AbstractNLSolver, rate_prototype)
    _du_cache = du_cache(nlsolver)
    return if _du_cache === nothing
        zero(rate_prototype)
    else
        first(_du_cache)
    end
end

mutable struct DAEResidualJacobianWrapper{
        isAD, F, pType, duType, uType, alphaType,
        gammaType,
        tmpType, uprevType, tType,
    } <: Function
    f::F
    p::pType
    tmp_du::duType
    tmp_u::uType
    α::alphaType
    invγdt::gammaType
    tmp::tmpType
    uprev::uprevType
    t::tType
    function DAEResidualJacobianWrapper(alg, f::F, p, α, invγdt, tmp, uprev, t) where {F}
        ad = ADTypes.dense_ad(alg_autodiff(alg))
        isautodiff = ad isa AutoForwardDiff
        if isautodiff
            tmp_du = DiffCache(uprev)
            tmp_u = DiffCache(uprev)
        else
            tmp_du = similar(uprev)
            tmp_u = similar(uprev)
        end
        return new{
            isautodiff, typeof(f), typeof(p), typeof(tmp_du), typeof(tmp_u), typeof(α),
            typeof(invγdt), typeof(tmp), typeof(uprev), typeof(t),
        }(
            f, p, tmp_du, tmp_u, α,
            invγdt, tmp, uprev, t
        )
    end
end

function ConstructionBase.setproperties(wrap::DAEResidualJacobianWrapper, patch::NamedTuple)
    for key in keys(patch)
        setproperty!(wrap, key, patch[key])
    end
    return wrap
end

is_autodiff(m::DAEResidualJacobianWrapper{isAD}) where {isAD} = isAD

function (m::DAEResidualJacobianWrapper)(out, x)
    if is_autodiff(m)
        tmp_du = get_tmp(m.tmp_du, x)
        tmp_u = get_tmp(m.tmp_u, x)
    else
        tmp_du = m.tmp_du
        tmp_u = m.tmp_u
    end
    @. tmp_du = (m.α * x + m.tmp) * m.invγdt
    @. tmp_u = x + m.uprev
    return m.f(out, tmp_du, tmp_u, m.p, m.t)
end

mutable struct DAEResidualDerivativeWrapper{
        F, pType, alphaType, gammaType, tmpType,
        uprevType, tType,
    }
    f::F
    p::pType
    α::alphaType
    invγdt::gammaType
    tmp::tmpType
    uprev::uprevType
    t::tType
end

function (m::DAEResidualDerivativeWrapper)(x)
    tmp_du = (m.α * x + m.tmp) * m.invγdt
    tmp_u = x + m.uprev
    return m.f(tmp_du, tmp_u, m.p, m.t)
end

SciMLBase.has_jac(f::DAEResidualJacobianWrapper) = SciMLBase.has_jac(f.f)
SciMLBase.has_Wfact(f::DAEResidualJacobianWrapper) = SciMLBase.has_Wfact(f.f)
SciMLBase.has_Wfact_t(f::DAEResidualJacobianWrapper) = SciMLBase.has_Wfact_t(f.f)

SciMLBase.has_jac(f::DAEResidualDerivativeWrapper) = SciMLBase.has_jac(f.f)
SciMLBase.has_Wfact(f::DAEResidualDerivativeWrapper) = SciMLBase.has_Wfact(f.f)
SciMLBase.has_Wfact_t(f::DAEResidualDerivativeWrapper) = SciMLBase.has_Wfact_t(f.f)

# Wrappers for computing separated DAE Jacobians: dF/du and dF/d(du).
# These pass the differentiation variable directly to f, so ForwardDiff
# handles mixed Float64/Dual arithmetic natively. No DiffCache needed.

# IIP: varies u only, du held fixed → Jacobian is dF/du
mutable struct DAEJacobianUWrapper{F, pType, duType, tType} <: Function
    f::F
    p::pType
    du_fixed::duType
    t::tType
end

function (m::DAEJacobianUWrapper)(out, u)
    return m.f(out, m.du_fixed, u, m.p, m.t)
end

# IIP: varies du only, u held fixed → Jacobian is dF/d(du)
mutable struct DAEJacobianDuWrapper{F, pType, uType, tType} <: Function
    f::F
    p::pType
    u_fixed::uType
    t::tType
end

function (m::DAEJacobianDuWrapper)(out, du)
    return m.f(out, du, m.u_fixed, m.p, m.t)
end

# OOP: varies u only → Jacobian is dF/du
mutable struct DAEJacobianUDerivativeWrapper{F, pType, duType, tType}
    f::F
    p::pType
    du_fixed::duType
    t::tType
end

(m::DAEJacobianUDerivativeWrapper)(u) = m.f(m.du_fixed, u, m.p, m.t)

# OOP: varies du only → Jacobian is dF/d(du)
mutable struct DAEJacobianDuDerivativeWrapper{F, pType, uType, tType}
    f::F
    p::pType
    u_fixed::uType
    t::tType
end

(m::DAEJacobianDuDerivativeWrapper)(du) = m.f(du, m.u_fixed, m.p, m.t)

SciMLBase.has_jac(::DAEJacobianUWrapper) = false
SciMLBase.has_jac(::DAEJacobianDuWrapper) = false
SciMLBase.has_jac(::DAEJacobianUDerivativeWrapper) = false
SciMLBase.has_jac(::DAEJacobianDuDerivativeWrapper) = false
SciMLBase.has_Wfact(::DAEJacobianUWrapper) = false
SciMLBase.has_Wfact(::DAEJacobianDuWrapper) = false
SciMLBase.has_Wfact_t(::DAEJacobianUWrapper) = false
SciMLBase.has_Wfact_t(::DAEJacobianDuWrapper) = false

# Cache structs bundling DAE-specific separated Jacobian state
struct DAEJacobiansCache{JduType, ufuType, ufduType, jcuType, jcduType}
    J_du::JduType
    uf_u::ufuType
    uf_du::ufduType
    jac_config_u::jcuType
    jac_config_du::jcduType
end

struct DAEJacobiansConstantCache{JduType, ufuType, ufduType}
    J_du::JduType
    uf_u::ufuType
    uf_du::ufduType
end

function _nlalg_with_linsolve(inner_alg, linsolve)
    linsolve === nothing && return inner_alg
    !hasfield(typeof(inner_alg), :descent) && return inner_alg
    descent = inner_alg.descent
    !hasfield(typeof(descent), :linsolve) && return inner_alg
    descent.linsolve !== nothing && return inner_alg
    return remake(inner_alg; descent = remake(descent; linsolve = linsolve))
end

"""
    build_nlsolver(alg, [nlalg,] u, uprev, p, t, dt, f, rate_prototype,
                   uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c, [α,]
                   iip, verbose) -> AbstractNLSolver

Construct the nonlinear solver object (an [`AbstractNLSolver`](@ref)) that an
implicit algorithm `alg` uses to solve its stage equations. `γ` and `c` are the
stage's diagonal coefficient and abscissa, `α` an optional scaling, `iip` the
in-place flag; the nonlinear-solver algorithm defaults to `alg.nlsolve`. Allocates
the appropriate cache (Newton `W`/factorization, functional/Anderson buffers, …)
for the chosen `nlalg`.
"""
function build_nlsolver(
        alg, u, uprev, p, t, dt, f::F, rate_prototype,
        ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, γ, c,
        iip, verbose
    ) where {F, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits,
        tTypeNoUnits, γ, c, 1, iip, verbose
    )
end

function build_nlsolver(
        alg, u, uprev, p, t, dt, f::F, rate_prototype,
        ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, γ, c, α,
        iip, verbose
    ) where {F, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return build_nlsolver(
        alg, alg.nlsolve, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, iip, verbose
    )
end

function daenlf(ztmp, z, p)
    tmp, ustep, γ, α, tstep, k, invγdt, _p, dt, f = p
    return _compute_rhs!(tmp, ztmp, ustep, γ, α, tstep, k, invγdt, _p, dt, f, z)[1]
end

function odenlf(ztmp, z, p)
    tmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f = p
    return _compute_rhs!(
        tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f, z
    )[1]
end

function build_nlsolver(
        alg,
        nlalg::Union{
            NLFunctional, NLAnderson, NLNewton, NonlinearSolveAlg,
            HomotopyNonlinearSolveAlg,
        },
        u, uprev, p, t, dt,
        f::F, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        γ, c, α,
        ::Val{true}, verbose
    ) where {
        F, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits,
    }
    #TODO
    #nlalg = SciMLBase.handle_defaults(alg, nlalg)
    # define unitless type
    uTolType = real(uBottomEltypeNoUnits)
    isdae = alg isa DAEAlgorithm

    # define fields of non-linear solver
    z = zero(u)
    tmp = zero(u)
    ztmp = zero(u)

    # build cache of non-linear solver
    ustep = zero(u)
    tstep = zero(t)
    k = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    atmp .= false
    dz = zero(u)

    if nlalg isa HomotopyNonlinearSolveAlg
        isdae && throw(
            ArgumentError(
                "HomotopyNonlinearSolveAlg does not support DAE problems: the step-size " *
                    "embedding degenerates the algebraic equations at λ = 0."
            )
        )
        invγdt = inv(oneunit(t) * one(uTolType))
        nlfunc = NonlinearFunction{true, SciMLBase.FullSpecialize}(homotopy_odenlf)
        nlcache = HomotopyNonlinearSolveCache(ustep, tstep, k, invγdt, nlfunc, Ref(0))
    elseif nlalg isa Union{NLNewton, NonlinearSolveAlg}
        nf = nlsolve_f(f, alg)

        # TODO: check if the solver is iterative
        weight = zero(u)

        if islinear(f) || SciMLBase.has_jac(f) ||
                (SciMLBase.has_jac_u(f) && SciMLBase.has_jac_du(f))
            du1 = rate_prototype
            uf = nothing
            jac_config = nothing
        else
            du1 = zero(rate_prototype)
            if isdae
                uf = DAEResidualJacobianWrapper(alg, f, p, α, inv(γ * dt), tmp, uprev, t)
            else
                uf = build_uf(alg, nf, t, p, Val(true))
            end
            jac_config = build_jac_config(alg, nf, uf, du1, uprev, u, ztmp, dz)
        end
        J, W = build_J_W(alg, u, uprev, p, t, dt, f, jac_config, uEltypeNoUnits, Val(true))

        tType = typeof(t)
        invγdt = inv(oneunit(t) * one(uTolType))

        if nlalg isa NonlinearSolveAlg
            γ = tTypeNoUnits(γ)
            α = tTypeNoUnits(α)
            dt = tTypeNoUnits(dt)
            W_for_reuse = if W isa WOperator && W.J !== nothing &&
                    !(W.J isa AbstractSciMLOperator)
                W._concrete_form
            else
                W
            end
            use_w_reuse = !isdae && f.nlstep_data === nothing &&
                W_for_reuse isa AbstractMatrix &&
                !(W_for_reuse isa AbstractSciMLOperator)
            prob = if f.nlstep_data !== nothing
                f.nlstep_data.nlprob
            else
                nlf = isdae ? daenlf : odenlf
                nlp_params = if isdae
                    (tmp, ustep, γ, α, tstep, k, invγdt, p, dt, f)
                else
                    (tmp, ustep, γ, α, tstep, k, invγdt, DIRK, p, dt, f)
                end
                if use_w_reuse
                    nlf_jac! = let W = W_for_reuse
                        (J_out, z, p) -> (copyto!(J_out, W); J_out)
                    end
                    # Without a `jac_prototype`, NonlinearSolve allocates a dense J for an
                    # analytic-jac function, so a sparse/structured W degrades to dense LU
                    # and sparse-only linsolves (e.g. KLU) crash. Mirror W's structure.
                    NonlinearProblem(
                        NonlinearFunction{true, SciMLBase.FullSpecialize}(
                            nlf;
                            jac = nlf_jac!,
                            jac_prototype = similar(W_for_reuse)
                        ),
                        ztmp, nlp_params
                    )
                else
                    NonlinearProblem(
                        NonlinearFunction{true, SciMLBase.FullSpecialize}(nlf),
                        ztmp, nlp_params
                    )
                end
            end
            inner_alg = _nlalg_with_linsolve(nlalg.alg, wrapprecs(alg.linsolve, W, weight))
            cache = init(prob, inner_alg, verbose = verbose.nonlinear_verbosity)
            nlcache = NonlinearSolveCache(
                ustep, tstep, k, atmp, invγdt, prob, cache,
                use_w_reuse ? J : nothing,
                use_w_reuse ? W_for_reuse : nothing,
                use_w_reuse ? uf : nothing,
                use_w_reuse ? jac_config : nothing,
                (use_w_reuse && uf !== nothing) ? du1 : nothing,
                weight,
                zero(tstep), true
            )
        else
            du = isdae ? k : nothing # k will be overwritten at solve time, but has the right type.
            linprob = LinearProblem(W, _vec(k), (du, u, p, t); u0 = _vec(dz))
            linsolve = init(
                linprob, wrapprecs(alg.linsolve, W, weight),
                alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
                assumptions = LinearSolve.OperatorAssumptions(true),
                verbose = verbose.linear_verbosity
            )

            # Build separated DAE Jacobian cache if applicable
            if isdae
                if islinear(f) || SciMLBase.has_jac(f) ||
                        (SciMLBase.has_jac_u(f) && SciMLBase.has_jac_du(f))
                    # User provides Jacobian(s) — no wrappers, but still need J_du storage
                    J_du = fill!(similar(J), false)
                    dae_jacobians = DAEJacobiansCache(
                        J_du, nothing, nothing, nothing, nothing
                    )
                else
                    uf_u = DAEJacobianUWrapper(f, p, zero(uprev), t)
                    uf_du = DAEJacobianDuWrapper(f, p, copy(uprev), t)
                    jac_config_u = build_jac_config(
                        alg, nf, uf_u, du1, uprev, u, ztmp, dz
                    )
                    jac_config_du = build_jac_config(
                        alg, nf, uf_du, du1, uprev, u, ztmp, dz
                    )
                    J_du = fill!(similar(J), false)
                    dae_jacobians = DAEJacobiansCache(
                        J_du, uf_u, uf_du, jac_config_u, jac_config_du
                    )
                end
            else
                dae_jacobians = nothing
            end

            nlcache = NLNewtonCache(
                ustep, tstep, k, atmp, dz, J, W, true,
                true, true, tType(dt), du1, uf, jac_config,
                linsolve, weight, invγdt, tType(nlalg.new_W_dt_cutoff), t,
                dae_jacobians
            )
        end
    elseif nlalg isa NLFunctional
        nlcache = NLFunctionalCache(ustep, tstep, k, atmp, dz)
    elseif nlalg isa NLAnderson
        max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
        Δz₊s = [zero(z) for i in 1:max_history]
        Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
        R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
        γs = Vector{uEltypeNoUnits}(undef, max_history)

        dzold = zero(z)
        z₊old = zero(z)

        nlcache = NLAndersonCache(
            ustep, tstep, atmp, k, dz, dzold, z₊old, Δz₊s, Q, R, γs,
            0,
            nlalg.aa_start, nlalg.droptol
        )
    end

    # build non-linear solver
    ηold = one(t)

    return NLSolver{true, tTypeNoUnits}(
        z, tmp, ztmp, tTypeNoUnits(γ), c, α, nlalg, nlalg.κ,
        nlalg.fast_convergence_cutoff, ηold, 0, nlalg.max_iter,
        Divergence, nlcache
    )
end

function oopdaenlf(z, p)
    tmp, α, tstep, invγdt, _p, dt, uprev, f = p
    return _compute_rhs(tmp, α, tstep, invγdt, p, dt, uprev, f, z)[1]
end

function oopodenlf(z, p)
    tmp, γ, α, tstep, invγdt, method, _p, dt, f = p
    return _compute_rhs(tmp, γ, α, tstep, invγdt, method, _p, dt, f, z)[1]
end

function build_nlsolver(
        alg,
        nlalg::Union{
            NLFunctional, NLAnderson, NLNewton, NonlinearSolveAlg,
            HomotopyNonlinearSolveAlg,
        },
        u, uprev, p,
        t, dt,
        f::F, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        γ, c, α,
        ::Val{false}, verbose
    ) where {
        F, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits,
    }
    #TODO
    #nlalg = SciMLBase.handle_defaults(alg, nlalg)
    # define unitless type
    uTolType = real(uBottomEltypeNoUnits)
    isdae = alg isa DAEAlgorithm

    # define fields of non-linear solver
    z = u
    tmp = u
    ztmp = u

    # build cache of non-linear solver
    tstep = zero(t)

    if nlalg isa HomotopyNonlinearSolveAlg
        isdae && throw(
            ArgumentError(
                "HomotopyNonlinearSolveAlg does not support DAE problems: the step-size " *
                    "embedding degenerates the algebraic equations at λ = 0."
            )
        )
        invγdt = inv(oneunit(t) * one(uTolType))
        nlfunc = NonlinearFunction{false, SciMLBase.FullSpecialize}(homotopy_oopodenlf)
        nlcache = HomotopyNonlinearSolveCache(nothing, tstep, nothing, invγdt, nlfunc, Ref(0))
    elseif nlalg isa Union{NLNewton, NonlinearSolveAlg}
        nf = nlsolve_f(f, alg)
        if isdae
            uf = DAEResidualDerivativeWrapper(f, p, α, inv(γ * dt), tmp, uprev, t)
        else
            uf = build_uf(alg, nf, t, p, Val(false))
        end

        tType = typeof(t)
        invγdt = inv(oneunit(t) * one(uTolType))

        J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(false))
        if nlalg isa NonlinearSolveAlg
            γ = tTypeNoUnits(γ)
            α = tTypeNoUnits(α)
            dt = tTypeNoUnits(dt)
            nlf = isdae ? oopdaenlf : oopodenlf
            nlp_params = if isdae
                (tmp, α, tstep, invγdt, p, dt, uprev, f)
            else
                (tmp, γ, α, tstep, invγdt, DIRK, p, dt, f)
            end
            prob = NonlinearProblem(
                NonlinearFunction{false, SciMLBase.FullSpecialize}(nlf),
                copy(ztmp), nlp_params
            )
            inner_alg = _nlalg_with_linsolve(nlalg.alg, alg.linsolve)
            cache = init(prob, inner_alg, verbose = verbose.nonlinear_verbosity)
            nlcache = NonlinearSolveCache(
                nothing, tstep, nothing, nothing, invγdt, prob, cache,
                nothing, nothing, nothing, nothing, nothing, nothing,
                zero(tstep), true
            )
        else
            # Build separated DAE Jacobian cache if applicable
            if isdae
                if SciMLBase.has_jac(f) ||
                        (SciMLBase.has_jac_u(f) && SciMLBase.has_jac_du(f))
                    dae_jacobians = DAEJacobiansConstantCache(
                        copy(J), nothing, nothing
                    )
                else
                    uf_u = DAEJacobianUDerivativeWrapper(f, p, zero(uprev), t)
                    uf_du = DAEJacobianDuDerivativeWrapper(f, p, copy(uprev), t)
                    dae_jacobians = DAEJacobiansConstantCache(
                        copy(J), uf_u, uf_du
                    )
                end
            else
                dae_jacobians = nothing
            end

            nlcache = NLNewtonConstantCache(
                tstep, J, W, true, true, true, tType(dt), uf,
                invγdt, tType(nlalg.new_W_dt_cutoff), t,
                dae_jacobians
            )
        end
    elseif nlalg isa NLFunctional
        nlcache = NLFunctionalConstantCache(tstep)
    elseif nlalg isa NLAnderson
        max_history = min(nlalg.max_history, nlalg.max_iter, length(z))
        Δz₊s = Vector{typeof(z)}(undef, max_history)
        Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
        R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
        γs = Vector{uEltypeNoUnits}(undef, max_history)

        dz = u
        dzold = u
        z₊old = u

        nlcache = NLAndersonConstantCache(
            tstep, dz, dzold, z₊old, Δz₊s, Q, R, γs, 0,
            nlalg.aa_start, nlalg.droptol
        )
    end

    # build non-linear solver
    ηold = one(tTypeNoUnits)
    return NLSolver{false, tTypeNoUnits}(
        z, tmp, ztmp, tTypeNoUnits(γ), c, α, nlalg, nlalg.κ,
        nlalg.fast_convergence_cutoff, ηold, 0, nlalg.max_iter,
        Divergence,
        nlcache
    )
end

## Anderson acceleration

"""
    anderson(z, cache)

Return the next iterate of the fixed-point iteration `z = g(z)` by performing Anderson
acceleration based on the current iterate `z` and the settings and history in the `cache`.
"""
@muladd function anderson(z, cache)
    (; dz, Δz₊s, z₊old, dzold, R, Q, γs, history, droptol) = cache

    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δz₊s)
    if history > max_history
        # circularly shift differences of z₊
        for i in 1:(max_history - 1)
            Δz₊s[i] = Δz₊s[i + 1]
        end

        # delete left-most column of QR decomposition
        qrdelete!(Q, R, max_history)

        # update size of history
        history = max_history
    end

    # update history of differences of z₊
    Δz₊s[history] = @.. broadcast = false z - z₊old

    # replace/add difference of residuals as right-most column to QR decomposition
    qradd!(Q, R, _vec(dz .- dzold), history)

    # update cached values
    cache.dzold = dz
    cache.z₊old = z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
        while cond(R) > droptol && history > 1
            qrdelete!(Q, R, history)
            history -= 1
            Qcur,
                Rcur = view(Q, :, 1:history),
                UpperTriangular(view(R, 1:history, 1:history))
        end
    end

    # save updated history
    cache.history = history

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))

    # update next iterate
    for i in 1:history
        z = @.. broadcast = false z - γs[i] * Δz₊s[i]
    end

    z
end

"""
    anderson!(z, cache)

Update the current iterate `z` of the fixed-point iteration `z = g(z)` in-place
by performing Anderson acceleration based on the settings and history in the `cache`.
"""
@muladd function anderson!(z, cache)
    (; dz, z₊old, dzold, Δz₊s, γs, R, Q, history, droptol) = cache

    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δz₊s)
    if history > max_history
        # circularly shift differences of z₊
        ptr = Δz₊s[1]
        for i in 1:(max_history - 1)
            Δz₊s[i] = Δz₊s[i + 1]
        end
        Δz₊s[max_history] = ptr

        # delete left-most column of QR decomposition
        qrdelete!(Q, R, max_history)

        # update size of history
        history = max_history
    end

    # update history of differences of z₊
    @.. broadcast = false Δz₊s[history] = z - z₊old

    # replace/add difference of residuals as right-most column to QR decomposition
    @.. broadcast = false dzold = dz - dzold
    qradd!(Q, R, _vec(dzold), history)

    # update cached values
    @.. broadcast = false dzold = dz
    @.. broadcast = false z₊old = z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
        while cond(R) > droptol && history > 1
            qrdelete!(Q, R, history)
            history -= 1
            Qcur,
                Rcur = view(Q, :, 1:history),
                UpperTriangular(view(R, 1:history, 1:history))
        end
    end

    # save updated history
    cache.history = history

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))

    # update next iterate
    for i in 1:history
        @.. broadcast = false z = z - γs[i] * Δz₊s[i]
    end

    nothing
end

## resize

function resize_nlsolver!(integrator::SciMLBase.DEIntegrator, i::Int)
    isdefined(integrator.cache, :nlsolver) || return

    (; nlsolver) = integrator.cache

    if nlsolver isa AbstractArray
        for idx in eachindex(nlsolver)
            resize!(nlsolver[idx], integrator, i)
        end
    else
        resize!(nlsolver, integrator, i)
    end

    nlsolver.alg isa NLNewton && resize!(nlsolver.cache.linsolve, i)

    # make it reset everything since the caches changed size!
    isnewton(nlsolver) && (nlsolver.cache.firstcall = true)

    return nothing
end

function Base.resize!(nlsolver::AbstractNLSolver, integrator, i::Int)
    resize!(nlsolver.z, i)
    resize!(nlsolver.tmp, i)
    resize!(nlsolver.ztmp, i)

    return resize!(nlsolver.cache, nlsolver, integrator, i)
end

## default: dispatch only on the cache
function Base.resize!(cache::AbstractNLSolverCache, nlsolver, integrator, i::Int)
    return Base.resize!(cache, i)
end

"""
    qrdelete!(Q, R, k)

Delete the left-most column of F = Q[:, 1:k] * R[1:k, 1:k] by updating Q and R.
Only Q[:, 1:(k-1)] and R[1:(k-1), 1:(k-1)] are valid on exit.
"""
function qrdelete!(Q::AbstractMatrix, R::AbstractMatrix, k::Int)
    n, m = size(Q)
    m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
    1 ≤ k ≤ m || throw(ArgumentError())

    # apply Givens rotations
    for i in 2:k
        g = first(givens(R, i - 1, i, i))
        lmul!(g, R)
        rmul!(Q, g')
    end

    # move columns of R
    @inbounds for j in 1:(k - 1)
        for i in 1:(k - 1)
            R[i, j] = R[i, j + 1]
        end
    end

    return Q, R
end

"""
    qradd!(Q, R, v, k)

Replace the right-most column of F = Q[:, 1:k] * R[1:k, 1:k] with v by updating Q and R.
This implementation modifies vector v as well. Only Q[:, 1:k] and R[1:k, 1:k] are valid on
exit.
"""
function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::AbstractVector, k::Int)
    n, m = size(Q)
    n == length(v) || throw(DimensionMismatch())
    m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
    1 ≤ k ≤ m || throw(ArgumentError())

    @inbounds for i in 1:(k - 1)
        q = view(Q, :, i)
        r = dot(q, v)

        R[i, k] = r
        axpy!(-r, q, v)
    end

    @inbounds begin
        d = norm(v)
        R[k, k] = d
        @.. broadcast = false @view(Q[:, k]) = v / d
    end

    return Q, R
end

function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::Number, k::Int)
    1 == LinearAlgebra.checksquare(Q) == LinearAlgebra.checksquare(R) ||
        throw(DimensionMismatch())
    k == 1 || throw(ArgumentError())

    R[1, 1] = abs(v)
    Q[1, 1] = one(v)

    return Q, R
end
