# Get mass_matrix from function, defaulting to I if field doesn't exist
# This allows nlsolver to work with DiscreteFunction which lacks mass_matrix
get_mass_matrix(f) = hasproperty(f, :mass_matrix) ? f.mass_matrix : I

# Not every `AbstractSciMLFunction` has an `nlstep_data` field: a `DynamicalODEFunction`
# has none, so reading it outright threw before the first step. Mirrors the field check
# `get_mass_matrix` already uses for `DiscreteFunction`; it constant-folds on a concrete
# `f`, so the branch selecting `prob` stays resolved at compile time.
get_nlstep_data(f) = hasfield(typeof(f), :nlstep_data) ? f.nlstep_data : nothing

# Map a flat linear-solve result back onto the state container. Delegates to
# ArrayInterface.restructure (which preserves ArrayPartition / other wrappers that
# plain reshape collapses); Number states need a convert because restructure
# relies on `similar`.
@inline _restructure_state(template, x) = ArrayInterface.restructure(template, x)
@inline _restructure_state(template::Number, x) = oftype(template, x)

get_status(nlsolver::AbstractNLSolver) = nlsolver.status
get_new_W_γdt_cutoff(nlsolver::AbstractNLSolver) = nlsolver.cache.new_W_γdt_cutoff
# handle FIRK
get_new_W_γdt_cutoff(alg::NewtonAlgorithm) = alg.new_W_γdt_cutoff

"""
    nlsolvefail(nlsolver::AbstractNLSolver) -> Bool
    nlsolvefail(status::NLStatus) -> Bool

Return whether a nonlinear solve failed.

# Arguments

  - `nlsolver`: a nonlinear solver whose current outcome is stored in its `status`
    property.
  - `status`: an [`OrdinaryDiffEqCore.NLStatus`](@ref) to classify directly.

# Returns

`true` for a non-positive status (`SlowConvergence`, `VerySlowConvergence`, or
`Divergence`) and `false` for `Convergence` or `FastConvergence`. Sibling solver
packages use this predicate after [`nlsolve!`](@ref) to decide whether to abandon
the current integrator step.
"""
nlsolvefail(nlsolver::AbstractNLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(status::NLStatus) = Int8(status) <= 0

relax(nlsolver::AbstractNLSolver) = relax(nlsolver.alg)
relax(alg::NLNewton) = alg.relax
relax(_) = 0 // 1

isnewton(nlsolver::AbstractNLSolver) = isnewton(nlsolver.cache)
isnewton(::AbstractNLSolverCache) = false
isnewton(::Union{NLNewtonCache, NLNewtonConstantCache}) = true

"""
    can_smooth_est(nlsolver::AbstractNLSolver) -> Bool

Return whether `nlsolver` can reuse its current `W` linear solve for an implicit
method's smoothed error estimate.

This capability query does not expose the nonlinear solver's cache or
factorization representation. A solver-author caller must use the ordinary
embedded estimate when this function returns `false`.
"""
can_smooth_est(nlsolver::AbstractNLSolver) = can_smooth_est(nlsolver.cache)
can_smooth_est(cache::AbstractNLSolverCache) = isnewton(cache)
can_smooth_est(cache::NonlinearSolveCache) = cache.linsolve !== nothing

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
    markfirststage!(nlsolver::AbstractNLSolver) -> Nothing

Mark `nlsolver` as being on the first implicit stage of the current integrator
step.

This mutates the solver's stage state when its algorithm tracks one and is a
no-op otherwise. Predictor and `W`-reuse logic may query the marker through
[`OrdinaryDiffEqCore.isfirststage`](@ref). The marker is cleared by
[`nlsolve!`](@ref)'s postamble after a solve.
"""
function markfirststage!(nlsolver::AbstractNLSolver)
    setfirststage!(nlsolver, true)
    return nothing
end

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
    du_alias_or_new(nlsolver::AbstractNLSolver, rate_prototype)

Return a derivative buffer compatible with `rate_prototype`.

# Arguments

  - `nlsolver`: the nonlinear solver that may already own a reusable derivative
    workspace.
  - `rate_prototype`: the integrator's derivative prototype and the template for
    a newly allocated buffer.

# Returns

The solver-owned derivative workspace when one is available; otherwise
`zero(rate_prototype)`. The returned object is mutable solver workspace when it
aliases an existing buffer, so callers may use it while constructing their
integrator cache but must not retain it beyond the nonlinear solver's lifetime.
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
    return remake(inner_alg; descent = remake(descent; linsolve))
end

# Analytic jacobian that copies the reused ODE `W` into the inner solver's own buffer each
# evaluation, so its factorization refreshes when the ODE rewrites `W`. A bare operator
# `jac_prototype` aliasing `W` would instead leave that factorization stale (its
# `update_coefficients!` is a no-op), silently converging the inner Newton on an out-of-date
# factorization.
struct WReuseJac{W} <: Function
    W::Base.RefValue{W}
end
(j::WReuseJac)(J_out, z, p) = (copyto!(J_out, j.W[]); J_out)

# The reused ODE `W` as the jacobian for the inner NonlinearSolve: a matrix-free `WOperator`
# passed through as an operator `jac_prototype` (applied via `mul!`, Krylov); a concrete W
# reused via an analytic `WReuseJac`. Used at build time and after `resize!` so the inner
# `NonlinearFunction`'s concrete type is preserved.
# Termination mode for the inner solve. The tolerances above are zero, so no criterion can
# ever fire — the integrator decides convergence. NonlinearSolve's default
# (`AbsNormSafeBestTerminationMode`) would nonetheless run its full per-iteration bookkeeping:
# best-objective tracking with a state copy, a circular objective trace, and — because it
# defaults to `max_stalled_steps = 32` — a state-difference buffer and a second norm. All of it
# is dead weight here, so ask for the plain absolute-norm mode instead.
_inner_termination() = NonlinearSolveBase.AbsNormTerminationMode(Base.Fix1(maximum, abs))

# Tolerances for the inner NonlinearSolve. `nlsolve!` drives that cache one `step!` at a time
# and decides convergence itself with the integrator's weighted `κ`/`η` test, exactly as it does
# for `NLNewton`, so the inner solver must never terminate on its own first: its default is an
# absolute residual bound unrelated to the integrator's tolerances, and the residual carries a
# `1/(γ·dt)` factor, so once `dt` grows it is met immediately — after which `step!` becomes a
# no-op, every later outer iteration sees a zero increment, and the outer test reads that as a
# perfect solve and accepts an unconverged stage.
#
# Zeroing `abstol`/`reltol` disables that criterion, but those same values are forwarded to the
# descent's linear solver, where zero is an unreachable target that costs an iterative (Krylov)
# solve its Newton-direction accuracy. `linsolve_kwargs` is splatted after them, so it restores
# a usable tolerance for the linear solve alone.
_inner_lintol(::Type{T}) where {T} = eps(real(one(T)))^(4 // 5)

"""
    set_inner_linear_reltol!(nlcache, integrator)

Point the inner solver's iterative linear solve at the integrator's relative tolerance,
the way `dolinsolve` does for `NLNewton`.

The `linsolve_kwargs` tolerance above is a fixed `eps^(4//5)`, so a Krylov descent solve
runs to ~1e-13 however loose the integrator is — at `reltol = 1e-3` that is ten orders
tighter than the step can use, and it is paid in Jacobian-vector products. Direct solvers
keep the fixed tolerance: `set_linear_reltol!` skips anything that needs a concrete `A`,
for which there is no tolerance to set.
"""
function set_inner_linear_reltol!(nlcache, integrator)
    lincache = get_linear_cache(nlcache)
    lincache === nothing && return nothing
    opts = integrator.opts
    # A fixed-step run has no `reltol` to inherit; `NLNewton` uses `eps` there.
    reltol = opts.adaptive ? opts.reltol : eps(real(eltype(get_u(nlcache))))
    set_linear_reltol!(lincache, reltol)
    return nothing
end

function reuse_jac_kwargs(W)
    Wr = W isa WOperator && W.J !== nothing && !(W.J isa AbstractSciMLOperator) ?
        W._concrete_form : W
    return Wr isa AbstractSciMLOperator ? (; jac_prototype = Wr) :
        (; jac = WReuseJac(Ref(Wr)), jac_prototype = (Z = similar(Wr); fill!(Z, 0); Z))
end

"""
    build_nlsolver(alg, [nlalg,] u, uprev, p, t, dt, f, rate_prototype,
                   uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c, [α,]
                   iip, verbose) -> AbstractNLSolver

Construct the nonlinear solver used by an implicit integrator algorithm.

# Arguments

  - `alg`: the implicit ODE, DAE, or stochastic integrator algorithm that owns
    the stage equation and differentiation configuration.
  - `nlalg`: optional nonlinear-solver algorithm. When omitted, `alg.nlsolve` is
    used.
  - `u`, `uprev`: current- and previous-state prototypes used to size solver
    workspaces.
  - `p`, `t`, `dt`: problem parameters, initial time, and initial step size.
  - `f`: the ODE or DAE function used by the implicit stage equation.
  - `rate_prototype`: derivative prototype used to size rate workspaces.
  - `uEltypeNoUnits`, `uBottomEltypeNoUnits`, `tTypeNoUnits`: unitless scalar
    types chosen by the integrator cache constructor.
  - `γ`, `c`: diagonal stage coefficient and stage abscissa.
  - `α`: optional stage scaling; defaults to `1`.
  - `iip`: `Val(true)` for an in-place problem and `Val(false)` for an
    out-of-place problem.
  - `verbose`: the integrator's differential-equation verbosity configuration.

# Returns

An [`OrdinaryDiffEqCore.AbstractNLSolver`](@ref). Its concrete type and cache are
implementation details; solver packages should retain the returned object and
operate on it through the documented nonlinear-solver interfaces. The factory
may retain aliases to the supplied prototypes as solver-owned workspace.

# Failure behavior

Unsupported nonlinear-solver algorithms fail by dispatch. A homotopy nonlinear
solver requested for a DAE throws `ArgumentError`; failures from Jacobian,
linear-solver, or inner nonlinear-solver initialization propagate to the caller.

# Solver-author usage

Implicit algorithm cache constructors call this factory once with their real
algorithm and problem prototypes, then pass the returned object to
[`nlsolve!`](@ref) for each stage. For example, a DIRK cache constructor uses the
form

```julia
nlsolver = build_nlsolver(
    alg, u, uprev, p, t, dt, f, rate_prototype,
    uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits,
    gamma, stage_abscissa, iip, verbose
)
```
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

# A DAE stage never reuses `W`, so the inner solver differentiates this residual itself and a
# dual-number backend calls it with a dual `z`; the cache's stage buffers hold the problem's
# own element type and cannot take the duals.
@inline function dae_stage_buffers(
        ustep::AbstractArray{T}, dustep::AbstractArray{T}, ::AbstractArray{T}
    ) where {T}
    return ustep, dustep
end
@inline function dae_stage_buffers(ustep, dustep, resid)
    R = eltype(resid)
    return similar(ustep, R), similar(dustep, R)
end

# The DAE kernel writes `du` at the stage into its second argument and the residual into its
# sixth, the reverse of the ODE kernel's single output: `dustep` is scratch, and the residual
# goes to the buffer the inner solver owns.
function daenlf(resid, z, p)
    tmp, ustep, α, tstep, dustep, invγdt, _p, uprev, f = p
    us, dus = dae_stage_buffers(ustep, dustep, resid)
    return _compute_rhs!(tmp, dus, us, α, tstep, resid, invγdt, _p, uprev, f, z)[1]
end

function odenlf(ztmp, z, p)
    tmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f = p
    return _compute_rhs!(
        tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f, z
    )[1]
end

# ---------------------------------------------------------------------------------------
# Stage-coordinate adapters for the `precondition`/`postcondition` solver options
#
# The unknown of the stage `NonlinearProblem` is the increment `z`, not the ODE state: the
# state at the stage is `ustep = tmp + γ⋅z` (`DIRK`) or `ustep = z`
# (`COEFFICIENT_MULTISTEP`), and the problem's parameter object is the internal tuple
# assembled by `build_nlsolver`/`initialize!` rather than the ODE's `p`. A user's corrector
# is written against the state and the ODE parameters, so it is conjugated with that affine
# map instead of being forwarded verbatim — `H` acting on `z` would clip an increment as
# though it were a concentration.

"""
    StageConditioner{iip, kind}(h, ubuf, uprevbuf, tmp, γ, method, p, pred)

Adapter presenting the stage state `ustep` and the ODE parameters to a `precondition`
(`kind === :pre`) or `postcondition` (`kind === :post`) supplied to
[`NonlinearSolveAlg`](@ref), and mapping a corrected state back to the stage increment.

The stage data defining that map is held here and refreshed by
[`refresh_stage_conditioners!`](@ref) rather than read from the parameter object the
corrector is handed: `reinit!` updates the inner cache's `p` but not its `prob`, and it is
`prob.p` that reaches a corrector, so those parameters are the ones the cache was *built*
with — a `tmp` and a `method` from before the first step.
"""
mutable struct StageConditioner{iip, kind, H, B, uType, γType, pType}
    const h::H
    const ubuf::B
    const uprevbuf::B
    tmp::uType
    γ::γType
    method::MethodType
    p::pType
    # The corrected stage predictor, and whether the next in-loop call is the first of this
    # stage. `reinit!` refreshes the inner cache's iterate but not the previous-iterate
    # buffer it hands a corrector, so on the first commit after a `reinit!` that buffer
    # still holds the *previous* stage's last `z` — read through the current stage's `tmp`
    # it is not a state at all. The predictor is what "the previous iterate" means there.
    pred::uType
    at_predictor::Bool
    # `init` evaluates the residual and corrects the initial guess as soon as the cache is
    # built, which for the stage problem happens before any stage exists: `tmp` is still a
    # zero buffer and the iterate is a placeholder. Staying inert until the first
    # `initialize!` keeps a user's `G`/`H` from being called on that placeholder state.
    active::Bool
end

function StageConditioner{iip, kind}(
        h::H, ubuf::B, uprevbuf::B, tmp::uType, γ::γType, method, p::pType, pred::uType
    ) where {iip, kind, H, B, uType, γType, pType}
    return StageConditioner{iip, kind, H, B, uType, γType, pType}(
        h, ubuf, uprevbuf, tmp, γ, method, p, pred, false, false
    )
end

# `precondition`: only the iterate and the parameters are remapped. The residual stays the
# stage residual — it is what the solver actually has, and it has no reading as a state.
function (c::StageConditioner{false, :pre})(fu, z, _)
    c.active || return fu
    return c.h(fu, compute_ustep(c.tmp, c.γ, z, c.method), c.p)
end

function (c::StageConditioner{true, :pre})(fu, z, _)
    c.active || return nothing
    c.h(fu, compute_ustep!(c.ubuf, c.tmp, c.γ, z, c.method), c.p)
    return nothing
end

# The inverse map is applied to the *change* the corrector made, `z + (H(ustep) - ustep)/γ`,
# rather than by recomputing `(H(ustep) - tmp)/γ`. The two agree in exact arithmetic, but
# only the first returns `z` bit-for-bit when the corrector leaves the state alone — and a
# limiter is inactive exactly when the iteration is converging, where a one-ulp perturbation
# of `z` is not harmless: it is enough to stop the stage residual from being exactly zero,
# which in turn stops the inner cache (deliberately given zero tolerances) from ever
# terminating, and `nlsolve!` reads a zero displacement from a non-terminated cache as a
# rejected trial step rather than as convergence.
function (c::StageConditioner{false, :post})(z, zprev, _, cache)
    c.active || return z
    zprev = take_previous_iterate!(c, zprev)
    c.method === COEFFICIENT_MULTISTEP && return c.h(z, zprev, c.p, cache)
    tmp, γ = c.tmp, c.γ
    ustep = compute_ustep(tmp, γ, z, c.method)
    uprevstep = compute_ustep(tmp, γ, zprev, c.method)
    return z + (c.h(ustep, uprevstep, c.p, cache) - ustep) / γ
end

function (c::StageConditioner{true, :post})(z, zprev, _, cache)
    c.active || return z
    zprev = take_previous_iterate!(c, zprev)
    if c.method === COEFFICIENT_MULTISTEP
        c.h(z, zprev, c.p, cache)
        return z
    end
    ustep, uprevstep, tmp, γ = c.ubuf, c.uprevbuf, c.tmp, c.γ
    @.. ustep = tmp + γ * z
    @.. uprevstep = tmp + γ * zprev
    c.h(ustep, uprevstep, c.p, cache)
    # `tmp + γ * z` is recomputed rather than kept in a third buffer: it is the same
    # deterministic expression as above, so an untouched `ustep` cancels exactly.
    @.. z = z + (ustep - (tmp + γ * z)) / γ
    return z
end

# The first in-loop commit of a stage follows the predictor, not the previous stage's last
# iterate (see the `pred` field).
function take_previous_iterate!(c::StageConditioner, zprev)
    c.at_predictor || return zprev
    c.at_predictor = false
    return c.pred
end

set_stage_predictor!(::Any, z) = nothing
function set_stage_predictor!(c::StageConditioner{true, :post}, z)
    copyto!(c.pred, z)
    c.at_predictor = true
    return nothing
end
function set_stage_predictor!(c::StageConditioner{false, :post}, z)
    c.pred = z
    c.at_predictor = true
    return nothing
end

# `nothing` stays `nothing` so that `init` never sees an inert conditioning keyword.
build_stage_conditioner(::Val, ::Val, ::Nothing, u, tmp, γ, p) = nothing

function build_stage_conditioner(::Val{iip}, ::Val{kind}, h, u, tmp, γ, p) where {iip, kind}
    ubuf, uprevbuf = iip ? (zero(u), zero(u)) : (nothing, nothing)
    pred = iip ? zero(u) : u
    return StageConditioner{iip, kind}(h, ubuf, uprevbuf, tmp, γ, DIRK, p, pred)
end

resize_conditioner!(::Any, ::Int) = nothing
function resize_conditioner!(c::StageConditioner{true}, i::Int)
    resize!(c.ubuf, i)
    resize!(c.uprevbuf, i)
    resize!(c.pred, i)
    return nothing
end

refresh_stage_conditioner!(::Any, tmp, γ, method, p) = nothing
function refresh_stage_conditioner!(c::StageConditioner, tmp, γ, method, p)
    c.tmp = tmp
    c.γ = γ
    c.method = method
    c.p = p
    c.at_predictor = false
    c.active = true
    return nothing
end

"""
    refresh_stage_conditioners!(nlcache, tmp, γ, method, p)

Point the stage conditioners at the data of the stage about to be solved. Called from
`initialize!`, which is the one place that knows all four: `method` in particular flips
between `DIRK` and `COEFFICIENT_MULTISTEP` within a single BDF solve.
"""
function refresh_stage_conditioners!(nlcache::NonlinearSolveCache, tmp, γ, method, p)
    refresh_stage_conditioner!(nlcache.precondition, tmp, γ, method, p)
    refresh_stage_conditioner!(nlcache.postcondition, tmp, γ, method, p)
    return nothing
end

"""
    apply_stage_predictor!!(nlsolver) -> z

Apply the `postcondition` corrector once to the stage predictor, the stage analogue of the
`H(u0, u0, p)` correction NonlinearSolve applies to an initial guess when a cache is built.
`reinit!` does not redo that correction, so without this the first residual of every stage
would be evaluated at an uncorrected iterate. The corrected predictor is also what the
stage's first in-loop correction sees as its previous iterate.
"""
function apply_stage_predictor!!(nlsolver)
    post = nlsolver.cache.postcondition
    post === nothing && return nlsolver.z
    z = nlsolver.z = post(nlsolver.z, nlsolver.z, nothing, nothing)
    set_stage_predictor!(post, z)
    return z
end

# Keywords for the stage `init`. Splatting an empty NamedTuple when neither option is in
# use keeps the untouched path byte-identical to before.
function conditioning_kwargs(nlcache::NonlinearSolveCache)
    return conditioning_kwargs(nlcache.precondition, nlcache.postcondition)
end
conditioning_kwargs(::Nothing, ::Nothing) = (;)
conditioning_kwargs(pre, ::Nothing) = (; precondition = pre)
conditioning_kwargs(::Nothing, post) = (; postcondition = post)
conditioning_kwargs(pre, post) = (; precondition = pre, postcondition = post)

# `odenlf`/`daenlf` write the stage state and `f`'s output through the preallocated
# `Float64` buffers carried in the problem's parameter tuple, and an `AutoSpecialize` `f` is
# a `FunctionWrapper` accepting only `Float64` and OrdinaryDiffEq's own one-chunk duals, so
# the in-place stage residual cannot be evaluated at `ForwardDiff.Dual`. Normally nothing
# tries: the reused `W` is handed to the inner solver as its Jacobian. A `precondition`
# removes that reuse — `W` is the Jacobian of the raw residual, not of the composed one — so
# the inner solver differentiates the residual itself, and a dual-based AD backend fails deep
# inside ForwardDiff. Refuse up front instead.
function check_precondition_differentiable(nlalg::NonlinearSolveAlg, iip::Bool)
    (iip && nlalg.precondition !== nothing) || return nothing
    autodiff = hasproperty(nlalg.alg, :autodiff) ? nlalg.alg.autodiff : nothing
    # `nothing` is "let NonlinearSolve pick", and it picks ForwardDiff for these problems.
    isduals = autodiff === nothing || autodiff isa AutoForwardDiff ||
        autodiff isa ADTypes.AutoPolyesterForwardDiff
    isduals || return nothing
    throw(
        ArgumentError(
            "`precondition` on an in-place problem needs an inner algorithm that does not \
            build its Jacobian with ForwardDiff: it disables reuse of the integrator's `W` \
            (which is the Jacobian of the raw stage residual, not of the preconditioned \
            one), and the in-place stage residual writes through preallocated `Float64` \
            buffers that dual numbers cannot be stored in. Pass an explicit non-dual \
            backend, e.g. `NonlinearSolveAlg(NewtonRaphson(autodiff = AutoFiniteDiff()); \
            precondition = ...)`, or state the problem out-of-place, where the residual \
            allocates and every backend works."
        )
    )
end

check_precondition_differentiable(nlalg, iip::Bool) = nothing

# The stage system's relation to the ODE state is only known for the residuals built here.
function check_conditioning_supported(nlalg::NonlinearSolveAlg, f, isdae::Bool)
    (nlalg.precondition === nothing && nlalg.postcondition === nothing) && return nothing
    opt = nlalg.precondition !== nothing ? "precondition" : "postcondition"
    isdae && throw(
        ArgumentError(
            "`$(opt)` is not supported for `DAEProblem`s: the stage residual is the \
            implicit DAE residual and its unknown is the state increment, so the option's \
            state-space contract does not apply."
        )
    )
    (hasproperty(f, :nlstep_data) && f.nlstep_data !== nothing) && throw(
        ArgumentError(
            "`$(opt)` is not supported when the ODE function carries `nlstep_data`: the \
            stage system is generated by ModelingToolkit and OrdinaryDiffEq has no map \
            from its unknowns back to the ODE state to state the option in."
        )
    )
    return nothing
end

check_conditioning_supported(nlalg, f, isdae::Bool) = nothing

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
        γ = tTypeNoUnits(γ)
        α = tTypeNoUnits(α)
        invγdt = inv(oneunit(t) * one(uTolType))
        nlfunc = NonlinearFunction{true, SciMLBase.FullSpecialize}(homotopy_odenlf)
        nf = Ref(0)
        nlp_params = (tmp, ustep, γ, α, tstep, k, invγdt, DIRK, p, dt, f, nf)
        prob = SciMLBase.HomotopyProblem{true}(
            nlfunc, ztmp, nlp_params; λspan = (zero(γ), one(γ))
        )
        cache_abstol = nlalg.abstol === nothing ? zero(nlalg.κ * one(uTolType)) :
            nlalg.abstol
        cache = _init_homotopy_nonlinear_cache(
            prob, nlalg, cache_abstol, verbose.nonlinear_verbosity
        )
        nlcache = HomotopyNonlinearSolveCache(
            ustep, tstep, k, invγdt, nlfunc, nf, verbose.nonlinear_verbosity, cache, false
        )
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
            check_conditioning_supported(nlalg, f, isdae)
            check_precondition_differentiable(nlalg, true)
            γ = tTypeNoUnits(γ)
            α = tTypeNoUnits(α)
            dt = tTypeNoUnits(dt)
            precondition = build_stage_conditioner(
                Val(true), Val(:pre), nlalg.precondition, u, tmp, float(γ), p
            )
            postcondition = build_stage_conditioner(
                Val(true), Val(:post), nlalg.postcondition, u, tmp, float(γ), p
            )
            # A matrix-free W (a `WOperator` wrapping a `JVPCache`, built for a Krylov
            # linear solver) is reused as an operator: NonlinearSolve applies it via
            # `mul!` rather than rebuilding a residual-derived AD JVP.
            matrixfree_W = W isa WOperator && W.J isa AbstractSciMLOperator
            W_for_reuse = if W isa WOperator && W.J !== nothing &&
                    !(W.J isa AbstractSciMLOperator)
                W._concrete_form
            else
                W
            end
            # `W` is the Jacobian of the *raw* stage residual, so handing it to the inner
            # solver as the Jacobian of a preconditioned one would be a lie; let the inner
            # solver differentiate the composition it actually solves.
            nlstep_data = get_nlstep_data(f)
            use_w_reuse = !isdae && nlstep_data === nothing &&
                precondition === nothing &&
                (
                (
                    W_for_reuse isa AbstractMatrix &&
                        !(W_for_reuse isa AbstractSciMLOperator)
                ) || matrixfree_W
            )
            prob = if nlstep_data !== nothing
                nlstep_data.nlprob
            else
                nlf = isdae ? daenlf : odenlf
                nlp_params = if isdae
                    (tmp, ustep, α, tstep, k, invγdt, p, uprev, f)
                else
                    (tmp, ustep, γ, α, tstep, k, invγdt, DIRK, p, dt, f)
                end
                if use_w_reuse
                    NonlinearProblem(
                        NonlinearFunction{true, SciMLBase.FullSpecialize}(
                            nlf; reuse_jac_kwargs(W)...
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
            inner_alg = _nlalg_with_linsolve(
                nlalg.alg,
                wrapprecs(default_krylov_warm_start(alg.linsolve), W, weight)
            )
            # The integrator owns the convergence decision, exactly as it does for `NLNewton`:
            # `nlsolve!` drives this cache one `step!` at a time and stops on its own weighted
            # `κ`/`η` test. Zero tolerances keep the inner solver's own criterion from ever
            # firing first — otherwise it terminates on an absolute residual bound unrelated to
            # the integrator's tolerances (and scaled by `1/(γ·dt)`, so it is trivially met once
            # `dt` grows), after which `step!` becomes a no-op, every later outer iteration sees
            # `ndz == 0`, and the outer test reads that as perfect convergence and accepts the
            # step with an unconverged stage.
            cache = init(
                prob, inner_alg; verbose = verbose.nonlinear_verbosity,
                abstol = zero(uTolType), reltol = zero(uTolType),
                termination_condition = _inner_termination(),
                linsolve_kwargs = (;
                    abstol = _inner_lintol(uTolType), reltol = _inner_lintol(uTolType),
                ),
                conditioning_kwargs(precondition, postcondition)...
            )
            if cache isa NonlinearSolveNoInitCache
                # An inner algorithm with no `__init` (every SimpleNonlinearSolve algorithm)
                # lands in this fallback cache, which only records the kwargs and splats them
                # into a complete `solve` per outer iteration — it cannot be driven one
                # `step!` at a time, so the integrator cannot own its convergence. The zeroed
                # tolerances above would then make every inner solve run to `MaxIters`
                # (except when the residual lands on exactly 0.0), so rebuild without them:
                # the inner solver terminates on its own default (nonzero) tolerances.
                cache = init(
                    prob, inner_alg; verbose = verbose.nonlinear_verbosity,
                    conditioning_kwargs(precondition, postcondition)...
                )
            end
            # Smoothed estimate `W \ tmp` reuses the inner solver's own W factorization (see
            # NonlinearSolveCache); `nothing` when it exposes none (native scalar/StaticArray
            # solve) falls back to the raw estimate.
            est_linsolve = use_w_reuse ? get_linear_cache(cache) : nothing
            nlcache = NonlinearSolveCache(
                ustep, tstep, k, atmp, invγdt, prob, cache,
                use_w_reuse ? J : nothing,
                use_w_reuse ? W_for_reuse : nothing,
                use_w_reuse ? uf : nothing,
                use_w_reuse ? jac_config : nothing,
                (use_w_reuse && uf !== nothing) ? du1 : nothing,
                weight,
                dz,
                est_linsolve,
                zero(tstep), true, false, precondition, postcondition
            )
        else
            du = isdae ? k : nothing # k will be overwritten at solve time, but has the right type.
            linprob = LinearProblem(W, _vec(k), (du, u, p, t); u0 = _vec(dz))
            linsolver = default_krylov_warm_start(alg.linsolve)
            precs = if isnothing(linsolver)
                Pl = LinearSolve.InvPreconditioner(_WeightPreconditioner(_vec(weight)))
                Pr = _WeightPreconditioner(_vec(weight))
                (; Pl, Pr)
            else
                (;)
            end
            linsolve = init(
                linprob, wrapprecs(linsolver, W, weight);
                precs...,
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
    tmp, α, tstep, invγdt, _p, uprev, f = p
    return _compute_rhs(tmp, α, tstep, invγdt, _p, uprev, f, z)[1]
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
        γ = tTypeNoUnits(γ)
        α = tTypeNoUnits(α)
        invγdt = inv(oneunit(t) * one(uTolType))
        nlfunc = NonlinearFunction{false, SciMLBase.FullSpecialize}(homotopy_oopodenlf)
        nf = Ref(0)
        nlp_params = (tmp, γ, α, tstep, invγdt, DIRK, p, dt, f, nf)
        prob = SciMLBase.HomotopyProblem{false}(
            nlfunc, zero(z), nlp_params; λspan = (zero(γ), one(γ))
        )
        cache_abstol = nlalg.abstol === nothing ? zero(nlalg.κ * one(uTolType)) :
            nlalg.abstol
        cache = _init_homotopy_nonlinear_cache(
            prob, nlalg, cache_abstol, verbose.nonlinear_verbosity
        )
        nlcache = HomotopyNonlinearSolveCache(
            nothing, tstep, nothing, invγdt, nlfunc, nf,
            verbose.nonlinear_verbosity, cache, false
        )
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
            check_conditioning_supported(nlalg, f, isdae)
            γ = tTypeNoUnits(γ)
            α = tTypeNoUnits(α)
            dt = tTypeNoUnits(dt)
            precondition = build_stage_conditioner(
                Val(false), Val(:pre), nlalg.precondition, u, tmp, float(γ), p
            )
            postcondition = build_stage_conditioner(
                Val(false), Val(:post), nlalg.postcondition, u, tmp, float(γ), p
            )
            nlf = isdae ? oopdaenlf : oopodenlf
            nlp_params = if isdae
                (tmp, α, tstep, invγdt, p, uprev, f)
            else
                (tmp, γ, α, tstep, invγdt, DIRK, p, dt, f)
            end
            inner_alg = _nlalg_with_linsolve(nlalg.alg, alg.linsolve)
            prob = NonlinearProblem(
                NonlinearFunction{false, SciMLBase.FullSpecialize}(nlf),
                copy(ztmp), nlp_params
            )
            # Zero tolerances: the integrator owns convergence (see the in-place branch above).
            cache = init(
                prob, inner_alg; verbose = verbose.nonlinear_verbosity,
                abstol = zero(uTolType), reltol = zero(uTolType),
                termination_condition = _inner_termination(),
                linsolve_kwargs = (;
                    abstol = _inner_lintol(uTolType), reltol = _inner_lintol(uTolType),
                ),
                conditioning_kwargs(precondition, postcondition)...
            )
            W_ref = nothing
            if cache isa NonlinearSolveNoInitCache
                if !isdae && f.nlstep_data === nothing && W isa StaticWOperator
                    W_ref = Ref(W.W)
                    nlf_jac = let W_ref = W_ref
                        (z, p) -> W_ref[]
                    end
                    prob = NonlinearProblem(
                        NonlinearFunction{false, SciMLBase.FullSpecialize}(
                            nlf; jac = nlf_jac, jac_prototype = W.W
                        ),
                        copy(ztmp), nlp_params
                    )
                end
                # `solve!`-driven fallback cache: it must keep terminating on its own
                # default (nonzero) tolerances (see the in-place branch above).
                cache = init(
                    prob, inner_alg; verbose = verbose.nonlinear_verbosity,
                    conditioning_kwargs(precondition, postcondition)...
                )
            end
            nlcache = NonlinearSolveCache(
                nothing, tstep, nothing, nothing, invγdt, prob, cache,
                nothing, W_ref, W_ref === nothing ? nothing : uf,
                nothing, nothing, nothing, nothing, nothing,
                zero(tstep), true, false, precondition, postcondition
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
    anderson(z, cache) -> accelerated_z

Return an Anderson-accelerated iterate for the fixed-point iteration `z = g(z)`.

`z` is the current iterate and `cache` is the Anderson workspace initialized by
the calling nonlinear solver. The function updates the workspace's iteration
history but returns the accelerated state rather than mutating `z`. Solver
packages that reuse this helper own the workspace's construction and lifetime;
no concrete Anderson cache type is part of this API.

Linear-algebra failures while updating or solving the history least-squares
system propagate to the caller.
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
    anderson!(z, cache) -> Nothing

Update the current iterate `z` of the fixed-point iteration `z = g(z)` in place
using Anderson acceleration.

`cache` is the Anderson workspace initialized by the calling nonlinear solver.
Both `z` and the workspace history are mutated. Solver packages that reuse this
helper own the workspace's construction and lifetime; no concrete Anderson cache
type is part of this API. Linear-algebra failures propagate to the caller.
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
