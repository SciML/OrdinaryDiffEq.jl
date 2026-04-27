abstract type AbstractController end

"""
    CommonControllerOptions{T}(; qmin=nothing, qmax=nothing,
                                qmax_first_step=nothing,
                                gamma=nothing,
                                qsteady_min=nothing, qsteady_max=nothing,
                                failfactor=nothing)

Composable holder for the standard step-size knobs every adaptive
controller uses:

- `qmin` / `qmax`: lower / upper bounds on the per-step shrink/grow factor.
- `qmax_first_step`: looser `qmax` applied to the very first accepted step
  (mirrors Sundials CVODE — the initial dt from automatic step-size
  selection is approximate, so a much larger growth is allowed once).
- `gamma`: safety factor applied to the controller's predicted dt change.
- `qsteady_min` / `qsteady_max`: deadband — if the proposed factor lies
  inside this interval, dt is held constant.
- `failfactor`: post-Newton-failure shrink factor used by
  [`post_newton_controller!`](@ref).

Concrete controllers (`IController`, `PIController`, `PIDController`,
`PredictiveController`, plus algorithm-specific controllers like
`BDFController` / `JVODEController`) embed a `CommonControllerOptions`
and forward the corresponding accessors (`get_qmin`, `get_qmax`,
`get_qmax_first_step`, `get_gamma`, `get_qsteady_min`, `get_qsteady_max`,
`get_failfactor`) to it.

The single `T` parameter is the eventual element type of the resolved
options (typically `Float64`); during construction each field is
`Union{Nothing, T}`, with `nothing` denoting "use the algorithm-specific
default". [`resolve_basic`](@ref) is called at `setup_controller_cache`
time to fill nothings via `qmin_default(alg)`, `qmax_default(alg)`, …
so a user-constructed `BDFController()` picks up BDF-tuned defaults
(`qmax = 5`, `qsteady_min = 9//10`, `qsteady_max = 12//10`, …) while a
user-constructed `IController()` falls back to the generic defaults
(`qmin = 1//5`, `qmax = 10`, …). This matches the historical behavior
of the per-algorithm step-size knobs that used to live on the
`OrdinaryDiffEq` algorithm structs themselves.
"""
struct CommonControllerOptions{T}
    qmin::Union{Nothing, T}
    qmax::Union{Nothing, T}
    qmax_first_step::Union{Nothing, T}
    gamma::Union{Nothing, T}
    qsteady_min::Union{Nothing, T}
    qsteady_max::Union{Nothing, T}
    failfactor::Union{Nothing, T}
    # Explicit inner constructor — suppresses the default
    # `CommonControllerOptions(args...) where {T}` outer that Julia would
    # otherwise generate, which dispatches to `T = Nothing` when all args
    # are `nothing` and breaks resolution downstream.
    function CommonControllerOptions{T}(
            qmin, qmax, qmax_first_step,
            gamma, qsteady_min, qsteady_max, failfactor,
        ) where {T}
        return new{T}(
            qmin, qmax, qmax_first_step,
            gamma, qsteady_min, qsteady_max, failfactor,
        )
    end
end

# T-inferring outer constructor for the all-positional form.
function CommonControllerOptions(
        qmin, qmax, qmax_first_step, gamma, qsteady_min, qsteady_max, failfactor,
    )
    T = _common_options_T(
        qmin, qmax, qmax_first_step,
        gamma, qsteady_min, qsteady_max, failfactor
    )
    return CommonControllerOptions{T}(
        qmin, qmax, qmax_first_step, gamma, qsteady_min, qsteady_max, failfactor,
    )
end

function CommonControllerOptions(;
        qmin = nothing, qmax = nothing, qmax_first_step = nothing,
        gamma = nothing, qsteady_min = nothing, qsteady_max = nothing,
        failfactor = nothing,
    )
    return CommonControllerOptions(
        qmin, qmax, qmax_first_step, gamma, qsteady_min, qsteady_max, failfactor,
    )
end

# Pick the element type T from the first concretely-set field; fall back
# to `Float64`. Step-size knobs need to hold rational defaults like
# `1//5`, so an integer (when the user passed a plain `Int` like
# `qmax = 3`) is promoted to a float type before resolving.
@inline _common_options_T() = Float64
@inline _common_options_T(::Nothing, rest...) = _common_options_T(rest...)
@inline _common_options_T(x, _...) = _floatify_QT(typeof(x))
@inline _floatify_QT(::Type{T}) where {T <: Integer} = float(T)
@inline _floatify_QT(::Type{T}) where {T} = T

"""
    resolve_basic(opts::CommonControllerOptions, alg, ::Type{QT})

Return a fully-resolved `CommonControllerOptions{QT}` where every
`nothing` field in `opts` has been replaced by the algorithm-specific
default (`qmin_default(alg)`, `qmax_default(alg)`, …). Called from each
controller's `setup_controller_cache` method so that user-constructed
controllers (e.g. `BDFController(qmax = 20)`) get sensible per-algorithm
defaults for the knobs they didn't set.
"""
function resolve_basic(opts::CommonControllerOptions, alg, ::Type{QT}) where {QT}
    @inline _resolve(::Nothing, default) = QT(default)
    @inline _resolve(value, _) = QT(value)
    return CommonControllerOptions{QT}(
        _resolve(opts.qmin, qmin_default(alg)),
        _resolve(opts.qmax, qmax_default(alg)),
        _resolve(opts.qmax_first_step, qmax_first_step_default(alg)),
        _resolve(opts.gamma, gamma_default(alg)),
        _resolve(opts.qsteady_min, qsteady_min_default(alg)),
        _resolve(opts.qsteady_max, qsteady_max_default(alg)),
        _resolve(opts.failfactor, failfactor_default(alg)),
    )
end

# Pick a sensible QT (element type) for resolving an unresolved
# `CommonControllerOptions` when none was provided explicitly.
@inline _resolved_QT(opts::CommonControllerOptions{T}) where {T} = _floatify_QT(T)

"""
    AbstractControllerCache

Each controller cache is expected to own the scalar error estimate used by the
step-size logic and exposed to users as `get_EEst(integrator)`. The accessors
[`get_EEst`](@ref) and [`set_EEst!`](@ref) read/write that scalar; the default
implementations dispatch on a `EEst` field. Controllers that track multiple
error estimates (e.g. BDF methods with `EEst1`, `EEst2`) can either still hold
a canonical scalar EEst or override the accessors directly.
"""
abstract type AbstractControllerCache end

"""
    get_EEst(controller_cache)
    get_EEst(integrator)

Return the scalar error estimate stored by the controller cache. When passed
an integrator, forwards to its `controller_cache`. The fallback reads the
`EEst` field on the cache; controllers that represent the error estimate
differently can override this method.
"""
@inline get_EEst(cache::AbstractControllerCache) = getfield(cache, :EEst)
@inline get_EEst(integrator::SciMLBase.DEIntegrator) =
    get_EEst(getfield(integrator, :controller_cache))

"""
    set_EEst!(controller_cache, val)
    set_EEst!(integrator, val)

Store `val` as the scalar error estimate on the controller cache. When passed
an integrator, forwards to its `controller_cache`. The fallback writes the
`EEst` field on the cache; controllers that represent the error estimate
differently can override this method.
"""
@inline set_EEst!(cache::AbstractControllerCache, val) =
    setfield!(cache, :EEst, convert(fieldtype(typeof(cache), :EEst), val))
@inline set_EEst!(integrator::SciMLBase.DEIntegrator, val) =
    set_EEst!(getfield(integrator, :controller_cache), val)

"""
    setup_controller_cache(alg, algcache, controller::AbstractController, ::Type{EEstT})::AbstractControllerCache

This function takes a controller together with the time stepping algorithm to
construct and initialize the respective cache for the controller. The
`EEstT` type parameter is the element type of the error estimate stored on
the returned cache (matches what used to live on `get_EEst(integrator)`).
"""
setup_controller_cache

# Back-compat dispatch for any 3-arg implementation still in the ecosystem.
# New code should pass the `EEstT` type parameter explicitly.
setup_controller_cache(alg, cache, controller::AbstractController) =
    setup_controller_cache(alg, cache, controller, Float64)

"""
    accept_step_controller(integrator, alg)::Bool
    accept_step_controller(integrator, controller_cache::AbstractControllerCache, alg)::Bool

This function decides whether the current time step should be accepted or rejected.
A return value of `false` corresponds to a rejection.
"""
accept_step_controller

"""
    stepsize_controller!(integrator, alg)
    stepsize_controller!(integrator, controller_cache::AbstractControllerCache, alg)

Update the cache to compute the new order of the time marching algorithm and prepare
for an update of the time step. The update can be either due to a rejection or acceptance
of the current time step.
"""
stepsize_controller!

"""
    step_accept_controller!(integrator, alg, q)
    step_accept_controller!(integrator, controller_cache::AbstractControllerCache, alg, q)

This function gets called in case of an accepted time step right after [`stepsize_controller!`](@ref).
It returns the proposed new time step length. Please note that the time step length might not be
applied as is and subject to further modification to e.g. match the next time stop.
"""
step_accept_controller!

"""
    step_reject_controller!(integrator, alg)
    step_reject_controller!(integrator, controller_cache::AbstractControllerCache, alg)

This function gets called in case of a rejected time step right after [`stepsize_controller!`](@ref).
It directly sets the time step length (i.e. `integrator.dt`).
"""
step_reject_controller!

# checks whether the controller should accept a step based on the error estimate
@inline function accept_step_controller(integrator, alg)
    return accept_step_controller(integrator, integrator.controller_cache, alg)
end
@inline function accept_step_controller(integrator, cache::Union{AbstractControllerCache, OrdinaryDiffEqCache}, alg)
    return get_EEst(integrator) <= 1
end

@inline function stepsize_controller!(integrator, alg)
    return stepsize_controller!(integrator, integrator.controller_cache, alg)
end
# Current fallback. This should actually dispatch onto the algorithms caches controller cache
@inline function stepsize_controller!(integrator, cache::OrdinaryDiffEqCache, alg)
    return stepsize_controller!(integrator, integrator.controller_cache, alg)
end

@inline function step_accept_controller!(integrator, alg, q)
    return step_accept_controller!(integrator, integrator.controller_cache, alg, q)
end

@inline function step_reject_controller!(integrator, alg)
    step_reject_controller!(integrator, integrator.controller_cache, alg)
    cache = integrator.cache
    if hasfield(typeof(cache), :nlsolve)
        nlsolve = cache.nlsolve
        nlsolve.prev_θ = one(nlsolve.prev_θ)
    end
    return nothing
end

"""
    get_current_qmax(integrator, qmax)

Return the effective maximum step size growth factor for the current step.
On the first step (before any successful steps), returns `qmax_first_step`
from the controller (defaulting to 10000). This allows a much larger step
size increase on the first step since the initial dt from the automatic
step size selection algorithm is only approximate.

This mirrors the behavior of Sundials CVODE, which limits h_new/h_old to 10
on normal steps but 10^4 on the first step.

See also: https://github.com/SciML/DifferentialEquations.jl/issues/299
"""
@inline function get_current_qmax(integrator, qmax)
    if integrator.success_iter == 0
        return get_qmax_first_step(integrator)
    end
    return qmax
end

"""
    get_qmin(integrator)
    get_qmax(integrator)
    get_qmax_first_step(integrator)
    get_gamma(integrator)
    get_qsteady_min(integrator)
    get_qsteady_max(integrator)
    get_failfactor(integrator)

Read a step-size knob from the integrator's controller. Default
dispatch reads `integrator.controller_cache.controller.basic.X` —
i.e. it goes through the `CommonControllerOptions` embedded on every concrete
controller (`IController`/`PIController`/`PIDController`/
`PredictiveController`/`BDFController`/`JVODEController`).

`CompositeControllerCache` overrides each accessor to delegate to the
currently active sub-cache (mirroring how `stepsize_controller!` and
friends dispatch). The transitional `DummyControllerCache` also
provides overrides for the BDF/Nordsieck cases that haven't been
migrated yet.

These accessors are what the integrator-level paths (e.g.
[`handle_step_rejection!`](@ref) for `qmin`,
[`post_newton_controller!`](@ref) for `failfactor`) call instead of
reading `integrator.opts.X` — the v7 controller refactor moved these
knobs off `DEOptions` and onto the controller object.
"""
@inline get_qmin(integrator::SciMLBase.DEIntegrator) =
    get_qmin(integrator, integrator.controller_cache)
@inline get_qmax(integrator::SciMLBase.DEIntegrator) =
    get_qmax(integrator, integrator.controller_cache)
@inline get_qmax_first_step(integrator::SciMLBase.DEIntegrator) =
    get_qmax_first_step(integrator, integrator.controller_cache)
@inline get_gamma(integrator::SciMLBase.DEIntegrator) =
    get_gamma(integrator, integrator.controller_cache)
@inline get_qsteady_min(integrator::SciMLBase.DEIntegrator) =
    get_qsteady_min(integrator, integrator.controller_cache)
@inline get_qsteady_max(integrator::SciMLBase.DEIntegrator) =
    get_qsteady_max(integrator, integrator.controller_cache)
@inline get_failfactor(integrator::SciMLBase.DEIntegrator) =
    get_failfactor(integrator, integrator.controller_cache)

# Default dispatch: reach through the controller's CommonControllerOptions.
@inline get_qmin(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.qmin
@inline get_qmax(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.qmax
@inline get_qmax_first_step(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.qmax_first_step
@inline get_gamma(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.gamma
@inline get_qsteady_min(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.qsteady_min
@inline get_qsteady_max(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.qsteady_max
@inline get_failfactor(integrator, cache::AbstractControllerCache) =
    cache.controller.basic.failfactor

reset_alg_dependent_opts!(controller::AbstractControllerCache, alg1, alg2) = nothing

reinit_controller!(integrator::SciMLBase.DEIntegrator, controller::AbstractControllerCache) = nothing
sync_controllers!(::AbstractControllerCache, ::AbstractControllerCache) = nothing

# Remember: Caches can also hold the control algorithm (see e.g. BDF and Nordsieck methods).
reinit_controller!(integrator::SciMLBase.DEIntegrator, cache::OrdinaryDiffEqCache) = nothing
sync_controllers!(::OrdinaryDiffEqCache, ::OrdinaryDiffEqCache) = nothing

function post_newton_controller!(integrator, alg)
    return post_newton_controller!(integrator, integrator.controller_cache, alg)
end
function post_newton_controller!(integrator, controller, alg)
    integrator.dt = integrator.dt / get_failfactor(integrator)
    return nothing
end

# This is a helper struct for algorithms with integrated controllers like Nordsieck and BDF methods.
struct DummyController <: AbstractController
end

"""
    DummyControllerCache

Controller cache used by algorithms that manage step-size selection themselves
(BDF, Nordsieck, Leaping, …). Holds the scalar error estimate exposed through
`get_EEst(integrator)` and a reference to the algorithm cache so existing dispatch
on the algorithm cache continues to work.
"""
mutable struct DummyControllerCache{T, C} <: AbstractControllerCache
    EEst::T
    cache::C
end

function setup_controller_cache(alg, cache, controller::DummyController, ::Type{E}) where {E}
    return DummyControllerCache{E, typeof(cache)}(oneunit(E), cache)
end

# Algorithms with integrated controllers (BDF, Nordsieck, …) only define their
# own `stepsize_controller!(integrator, alg)` 2-arg method that reaches for
# `integrator.cache`. When such an algorithm appears as a branch of a
# `CompositeAlgorithm`, the composite dispatch hands us the sub-cache
# (`DummyControllerCache`) explicitly, so fall back to the alg-level method.
@inline stepsize_controller!(integrator, ::DummyControllerCache, alg) =
    stepsize_controller!(integrator, alg)
@inline step_accept_controller!(integrator, ::DummyControllerCache, alg, q) =
    step_accept_controller!(integrator, alg, q)
@inline step_reject_controller!(integrator, ::DummyControllerCache, alg) =
    step_reject_controller!(integrator, alg)
@inline post_newton_controller!(integrator, ::DummyControllerCache, alg) =
    post_newton_controller!(integrator, alg)
@inline accept_step_controller(integrator, cache::DummyControllerCache, alg) =
    get_EEst(cache) <= 1
# DummyControllerCache is used by some SDE algorithms (e.g.
# StochasticDiffEqLeaping) that haven't been migrated to a proper
# CommonControllerOptions-based controller yet. They keep the step-size knobs as
# fields on the algorithm itself, so accessors fall back to those fields
# when present.
for (accessor, default) in (
        (:get_qmin, :(qmin_default(integrator.alg))),
        (:get_qmax, :(qmax_default(integrator.alg))),
        (:get_qmax_first_step, :(qmax_first_step_default(integrator.alg))),
        (:get_gamma, :(gamma_default(integrator.alg))),
        (:get_qsteady_min, :(qsteady_min_default(integrator.alg))),
        (:get_qsteady_max, :(qsteady_max_default(integrator.alg))),
        (:get_failfactor, :(failfactor_default(integrator.alg))),
    )
    field = Symbol(string(accessor)[5:end]) # strip leading "get_"
    @eval @inline function $accessor(integrator, ::DummyControllerCache)
        alg = integrator.alg
        return hasfield(typeof(alg), $(QuoteNode(field))) ?
            getfield(alg, $(QuoteNode(field))) : $default
    end
end


# Standard integral (I) step size controller
"""
    IController()

The standard (integral) controller is the most basic step size controller.
This controller is usually the first one introduced in numerical analysis classes
but should only be used rarely in practice because of efficiency problems for
many problems/algorithms.

Construct an integral (I) step size controller adapting the time step
based on the formula

```
Δtₙ₊₁ = εₙ₊₁^(1/k) * Δtₙ
```

where `k = get_current_adaptive_order(alg, integrator.cache) + 1` and `εᵢ` is the
inverse of the error estimate `get_EEst(integrator)` scaled by the tolerance
(Hairer, Nørsett, Wanner, 2008, Section II.4).
The step size factor is multiplied by the safety factor `gamma` and clipped to
the interval `[qmin, qmax]`.
A step will be accepted whenever the estimated error `get_EEst(integrator)` is
less than or equal to unity. Otherwise, the step is rejected and re-tried with
the predicted step size.

## References

  - Hairer, Nørsett, Wanner (2008)
    Solving Ordinary Differential Equations I Nonstiff Problems
    [DOI: 10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
"""
struct IController{B <: CommonControllerOptions} <: AbstractController
    basic::B
end

# Keyword form: holds nothings until `setup_controller_cache` resolves them
# against the algorithm.
IController(; kwargs...) = IController(CommonControllerOptions(; kwargs...))

# alg-aware forms: resolve immediately to a fully-typed controller.
IController(alg; kwargs...) = IController(Float64, alg; kwargs...)
IController(::Type{QT}, alg; kwargs...) where {QT} =
    IController(resolve_basic(CommonControllerOptions(; kwargs...), alg, QT))

mutable struct IControllerCache{T, E} <: AbstractControllerCache
    controller::IController{CommonControllerOptions{T}}
    dtreject::T
    EEst::E
end

function setup_controller_cache(alg, cache, controller::IController, ::Type{E}) where {E}
    QT = _resolved_QT(controller.basic)
    resolved = IController(resolve_basic(controller.basic, alg, QT))
    T = QT
    return IControllerCache{T, E}(resolved, T(1 // 10^4), oneunit(E))
end

@inline function stepsize_controller!(integrator, cache::IControllerCache, alg)
    (; qmin, qmax, gamma) = cache.controller.basic
    qmax = get_current_qmax(integrator, qmax)
    EEst = DiffEqBase.value(get_EEst(integrator))

    if iszero(EEst)
        q = inv(qmax)
    else
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = fastpower(EEst, expo) / gamma
        @fastmath q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        # TODO: Shouldn't this be in `step_accept_controller!` as for the PI controller?
        cache.dtreject = DiffEqBase.value(integrator.dt) / q
    end
    return q
end

# TODO change signature to remove the q input
function step_accept_controller!(integrator, cache::IControllerCache, alg, q)
    (; qsteady_min, qsteady_max) = cache.controller.basic

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, cache::IControllerCache, alg)
    return integrator.dt = cache.dtreject
end

reinit_controller!(integrator::SciMLBase.DEIntegrator, cache::IControllerCache) = nothing

function sync_controllers!(cache1::IControllerCache, cache2::IControllerCache)
    cache1.dtreject = cache2.dtreject
    return nothing
end

# PI step size controller
"""
    PIController(beta1, beta2)

The proportional-integral (PI) controller is a widespread step size controller
with improved stability properties compared to the [`IController`](@ref).
This controller is the default for most algorithms in OrdinaryDiffEq.jl.

Construct a PI step size controller adapting the time step based on the formula

```
Δtₙ₊₁ = εₙ₊₁^β₁ * εₙ^β₂ * Δtₙ
```

where `εᵢ` are inverses of the error estimates scaled by the tolerance
(Hairer, Nørsett, Wanner, 2010, Section IV.2).
The step size factor is multiplied by the safety factor `gamma` and clipped to
the interval `[qmin, qmax]`.
A step will be accepted whenever the estimated error `get_EEst(integrator)` is
less than or equal to unity. Otherwise, the step is rejected and re-tried with
the predicted step size.

!!! note

    The coefficients `beta1, beta2` are not scaled by the order of the method,
    in contrast to the [`PIDController`](@ref). For the `PIController`, this
    scaling by the order must be done when the controller is constructed.

## References

  - Hairer, Nørsett, Wanner (2010)
    Solving Ordinary Differential Equations II Stiff and Differential-Algebraic Problems
    [DOI: 10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
  - Hairer, Nørsett, Wanner (2008)
    Solving Ordinary Differential Equations I Nonstiff Problems
    [DOI: 10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
"""
mutable struct PIController{B <: CommonControllerOptions, T} <: AbstractController # TODO remove the mutable once AitkenNeville is fixed
    basic::B
    beta1::T
    beta2::T
    qoldinit::T
end

# Two-positional-arg form (beta1, beta2 explicit). Keyword-only knobs go
# into the embedded CommonControllerOptions and get resolved at setup time.
function PIController(beta1::Real, beta2::Real; qoldinit = 1 // 10^4, kwargs...)
    T = typeof(beta1)
    basic = CommonControllerOptions(; kwargs...)
    return PIController{typeof(basic), T}(basic, T(beta1), T(beta2), T(qoldinit))
end

# alg-aware forms: resolve CommonControllerOptions immediately, default beta1/beta2/qoldinit.
function PIController(alg; kwargs...)
    return PIController(Float64, alg; kwargs...)
end

function PIController(
        ::Type{QT}, alg;
        beta1 = nothing, beta2 = nothing, qoldinit = nothing,
        kwargs...
    ) where {QT}
    beta2 = beta2 === nothing ? beta2_default(alg) : beta2
    beta1 = beta1 === nothing ? beta1_default(alg, beta2) : beta1
    qoldinit = qoldinit === nothing ? 1 // 10^4 : qoldinit
    basic = resolve_basic(CommonControllerOptions(; kwargs...), alg, QT)
    return PIController{typeof(basic), QT}(basic, QT(beta1), QT(beta2), QT(qoldinit))
end

mutable struct PIControllerCache{T, E} <: AbstractControllerCache
    controller::PIController{CommonControllerOptions{T}, T}
    # Cached εₙ₊₁^β₁
    q11::T
    # Previous EEst
    errold::T
    EEst::E
end

function setup_controller_cache(alg, cache, controller::PIController, ::Type{E}) where {E}
    QT = _resolved_QT(controller.basic)
    basic = resolve_basic(controller.basic, alg, QT)
    resolved = PIController{typeof(basic), QT}(
        basic, QT(controller.beta1), QT(controller.beta2), QT(controller.qoldinit),
    )
    T = QT
    return PIControllerCache{T, E}(
        resolved, one(T), T(resolved.qoldinit), oneunit(E),
    )
end

@inline function stepsize_controller!(integrator, cache::PIControllerCache, alg)
    (; errold, controller) = cache
    (; qmin, qmax, gamma) = controller.basic
    qmax = get_current_qmax(integrator, qmax)
    (; beta1, beta2) = controller
    EEst = DiffEqBase.value(get_EEst(integrator))

    if iszero(EEst)
        q = inv(qmax)
    else
        q11 = fastpower(EEst, beta1)
        q = q11 / fastpower(errold, beta2)
        cache.q11 = q11
        @fastmath q = clamp(q / gamma, inv(qmax), inv(qmin))
    end
    return q
end

function step_accept_controller!(integrator, cache::PIControllerCache, alg, q)
    (; controller) = cache
    (; qsteady_min, qsteady_max) = controller.basic
    qoldinit = controller.qoldinit
    EEst = DiffEqBase.value(get_EEst(integrator))

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    cache.errold = max(EEst, qoldinit)
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, cache::PIControllerCache, alg)
    (; controller, q11) = cache
    (; qmin, gamma) = controller.basic
    return integrator.dt /= min(inv(qmin), q11 / gamma)
end

function reinit_controller!(integrator::SciMLBase.DEIntegrator, cache::PIControllerCache{T}) where {T}
    cache.q11 = one(T)
    return cache.errold = T(cache.controller.qoldinit)
end

function sync_controllers!(cache1::PIControllerCache, cache2::PIControllerCache)
    cache1.q11 = cache2.q11
    cache1.errold = cache2.errold
    return nothing
end

# PID step size controller
"""
    PIDController(beta1, beta2, beta3=zero(beta1);
                  limiter=default_dt_factor_limiter,
                  accept_safety=0.81)

The proportional-integral-derivative (PID) controller is a generalization of the
[`PIController`](@ref) and can have improved stability and efficiency properties.

Construct a PID step size controller adapting the time step based on the formula

```
Δtₙ₊₁ = εₙ₊₁^(β₁/k) * εₙ^(β₂/k) * εₙ₋₁^(β₃/ k) * Δtₙ
```

where `k = min(alg_order, alg_adaptive_order) + 1` and `εᵢ` are inverses of
the error estimates scaled by the tolerance (Söderlind, 2003).
The step size factor is limited by the `limiter` with default value

```
limiter(x) = one(x) + atan(x - one(x))
```

as proposed by Söderlind and Wang (2006). A step will be accepted whenever the
predicted step size change is bigger than `accept_safety`. Otherwise, the step
is rejected and re-tried with the predicted step size.

Some standard controller parameters suggested in the literature are

| Controller | `beta1` | `beta2` | `beta3` |
|:---------- | -------:| -------:|:-------:|
| basic      | `1.00`  | `0.00`  | `0`     |
| PI42       | `0.60`  | `-0.20` | `0`     |
| PI33       | `2//3`  | `-1//3` | `0`     |
| PI34       | `0.70`  | `-0.40` | `0`     |
| H211PI     | `1//6`  | `1//6`  | `0`     |
| H312PID    | `1//18` | `1//9`  | `1//18` |

!!! note

    In contrast to the [`PIController`](@ref), the coefficients `beta1, beta2, beta3`
    are scaled by the order of the method. Thus, standard controllers such as PI42
    can use the same coefficients `beta1, beta2, beta3` for different algorithms.

!!! note

    In contrast to other controllers, the `PIDController` does not use the keyword
    arguments `qmin, qmax` to limit the step size change or the safety factor `gamma`.
    These common keyword arguments are replaced by the `limiter` and `accept_safety`
    to guarantee a smooth behavior (Söderlind and Wang, 2006).
    Because of this, a `PIDController` behaves different from a [`PIController`](@ref),
    even if `beta1, beta2` are adapted accordingly and `iszero(beta3)`.

## References

  - Söderlind (2003)
    Digital Filters in Adaptive Time-Stepping
    [DOI: 10.1145/641876.641877](https://doi.org/10.1145/641876.641877)
  - Söderlind, Wang (2006)
    Adaptive time-stepping and computational stability
    [DOI: 10.1016/j.cam.2005.03.008](https://doi.org/10.1016/j.cam.2005.03.008)
  - Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)    # limiter of the dt factor (before clipping)
"""
struct PIDController{B <: CommonControllerOptions, QT, Limiter} <: AbstractController
    basic::B
    beta::NTuple{3, QT} # controller coefficients
    accept_safety::QT   # accept a step if the predicted change of the step size
    # is bigger than this parameter
    limiter::Limiter    # limiter of the dt factor (before clipping)
end

@inline default_dt_factor_limiter(x) = one(x) + atan(x - one(x))

# Two/three-positional form: betas explicit, CommonControllerOptions-related knobs
# can be passed as kwargs and stay nothing until setup_controller_cache.
function PIDController(
        beta1::Real, beta2::Real, beta3::Real = zero(beta1);
        accept_safety = 0.81, limiter = default_dt_factor_limiter, kwargs...
    )
    beta = map(float, promote(beta1, beta2, beta3))
    QT = typeof(beta[1])
    basic = CommonControllerOptions(; kwargs...)
    return PIDController{typeof(basic), QT, typeof(limiter)}(
        basic, beta, QT(accept_safety), limiter,
    )
end

function PIDController(alg; kwargs...)
    return PIDController(Float64, alg; kwargs...)
end

function PIDController(
        ::Type{QT}, alg;
        beta = nothing, accept_safety = 0.81,
        limiter = default_dt_factor_limiter,
        kwargs...
    ) where {QT}
    if beta === nothing
        beta2 = QT(beta2_default(alg))
        beta1 = QT(beta1_default(alg, beta2))
        beta3 = QT(zero(beta1))
    else
        beta1, beta2, beta3 = beta
    end
    beta = map(float, promote(beta1, beta2, beta3))
    basic = resolve_basic(CommonControllerOptions(; kwargs...), alg, QT)
    return PIDController{typeof(basic), QT, typeof(limiter)}(
        basic, beta, QT(accept_safety), limiter,
    )
end

function Base.show(io::IO, controller::PIDController)
    return print(
        io, "PIDController(beta=", controller.beta,
        ", accept_safety=", controller.accept_safety,
        ", limiter=", controller.limiter,
        ")"
    )
end

mutable struct PIDControllerCache{T, Limiter, E} <: AbstractControllerCache
    controller::PIDController{CommonControllerOptions{T}, T, Limiter}
    err::Vector{T} # history of the error estimates
    dt_factor::T
    EEst::E
end

function reinit_controller!(integrator::SciMLBase.DEIntegrator, cache::PIDControllerCache{T}) where {T}
    cache.err = ones(T, 3)
    cache.dt_factor = one(T)
    return nothing
end

function setup_controller_cache(alg, cache, controller::PIDController, ::Type{E}) where {E}
    QT = _resolved_QT(controller.basic)
    basic = resolve_basic(controller.basic, alg, QT)
    resolved = PIDController{typeof(basic), QT, typeof(controller.limiter)}(
        basic, map(QT, controller.beta), QT(controller.accept_safety), controller.limiter,
    )
    err = ones(QT, 3)
    return PIDControllerCache{QT, typeof(controller.limiter), E}(
        resolved, err, one(QT), oneunit(E),
    )
end

@inline function stepsize_controller!(integrator, cache::PIDControllerCache, alg)
    (; controller) = cache
    beta1, beta2, beta3 = controller.beta

    EEst = DiffEqBase.value(get_EEst(integrator))

    # If the error estimate is zero, we can increase the step size as much as
    # desired. This additional check fixes problems of the code below when the
    # error estimates become zero
    # -> err1, err2, err3 become Inf
    # -> err1^positive_number * err2^negative_number becomes NaN
    # -> dt becomes NaN
    #
    # `EEst_min` is smaller than PETSC_SMALL used in the equivalent logic in PETSc.
    # For example, `eps(Float64) ≈ 2.2e-16` but `PETSC_SMALL ≈ 1.0e-10` for `double`.
    EEst_min = eps(typeof(EEst))
    # The code below is a bit more robust than
    # ```
    # if iszero(EEst)
    #   EEst = eps(typeof(EEst))
    # end
    # ```
    EEst = max(EEst, EEst_min)

    cache.err[1] = inv(EEst)
    err1, err2, err3 = cache.err

    k = min(alg_order(alg), alg_adaptive_order(alg)) + 1
    dt_factor = err1^(beta1 / k) * err2^(beta2 / k) * err3^(beta3 / k)
    if isnan(dt_factor)
        @warn "unlimited dt_factor" dt_factor err1 err2 err3 beta1 beta2 beta3 k
    end
    cache.dt_factor = controller.limiter(dt_factor)

    # Note: No additional limiting of the form
    #   dt_factor = max(qmin, min(qmax, dt_factor))
    # is necessary since the `limiter` should take care of that. The default limiter
    # ensures
    #   0.21 ≈ limiter(0) <= dt_factor <= limiter(Inf) ≈ 2.57
    # See Söderlind, Wang (2006), Section 6.
    return cache.dt_factor
end

@inline function accept_step_controller(integrator, cache::PIDControllerCache, alg)
    return cache.dt_factor >= cache.controller.accept_safety
end

function step_accept_controller!(integrator, cache::PIDControllerCache, alg, dt_factor)
    (; controller) = cache
    (; qsteady_min, qsteady_max) = controller.basic

    if qsteady_min <= inv(dt_factor) <= qsteady_max
        dt_factor = one(dt_factor)
    end
    @inbounds begin
        cache.err[3] = cache.err[2]
        cache.err[2] = cache.err[1]
    end
    return integrator.dt * dt_factor # new dt
end

function step_reject_controller!(integrator, cache::PIDControllerCache, alg)
    return integrator.dt *= cache.dt_factor
end

function sync_controllers!(cache1::PIDControllerCache, cache2::PIDControllerCache)
    cache1.err = cache2.err
    cache1.dt_factor = cache2.dt_factor
    return nothing
end

# Gustafsson predictive step size controller
"""
    PredictiveController()

The Gustafsson acceleration algorithm accelerates changes so that way algorithms
can more swiftly change to handle quick transients. This algorithm is thus
well-suited for stiff solvers where this can be expected, and is the default
for algorithms like the (E)SDIRK methods.

```julia
(; qmin, qmax, gamma) = controller
qmax = get_current_qmax(integrator, qmax)
niters = integrator.cache.nlsolver.iter
fac = min(gamma,
    (1 + 2 * integrator.cache.nlsolver.maxiters) * gamma /
    (niters + 2 * integrator.cache.nlsolver.maxiters))
expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
qtmp = fastpower(get_EEst(integrator), expo) / fac
@fastmath q = max(inv(qmax), min(inv(qmin), qtmp))
cache.qold = q
q
```

where `controller` and `cache` are the controller and its cache stored in
`integrator.controller_cache`. `niters` is the number of Newton iterations
which was required in the most recent step of the algorithm. Note that
these values are used differently depending on acceptance and rejectance.
When the step is accepted, the following logic is applied:

```julia
(; dtacc, erracc) = cache
(; qmin, qmax, gamma, qsteady_min, qsteady_max) = controller
qmax = get_current_qmax(integrator, qmax)
if integrator.success_iter > 0
    expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
    qgus = (dtacc / integrator.dt) * fastpower((get_EEst(integrator)^2) / erracc, expo)
    qgus = max(inv(qmax), min(inv(qmin), qgus / gamma))
    qacc = max(q, qgus)
else
    qacc = q
end
if qsteady_min <= qacc <= qsteady_max
    qacc = one(qacc)
end
cache.dtacc = integrator.dt
cache.erracc = max(1e-2, get_EEst(integrator))
integrator.dt / qacc
```

When it rejects, it's the same as the [`IController`](@ref):

```julia
if integrator.success_iter == 0
    integrator.dt *= 0.1
else
    integrator.dt = integrator.dt / cache.qold
end
```
"""
struct PredictiveController{B <: CommonControllerOptions} <: AbstractController
    basic::B
end

PredictiveController(; kwargs...) = PredictiveController(CommonControllerOptions(; kwargs...))

PredictiveController(alg; kwargs...) = PredictiveController(Float64, alg; kwargs...)
PredictiveController(::Type{QT}, alg; kwargs...) where {QT} =
    PredictiveController(resolve_basic(CommonControllerOptions(; kwargs...), alg, QT))

mutable struct PredictiveControllerCache{T, E} <: AbstractControllerCache
    controller::PredictiveController{CommonControllerOptions{T}}
    dtacc::T
    erracc::T
    qold::T
    EEst::E
end

function reinit_controller!(integrator::SciMLBase.DEIntegrator, cache::PredictiveControllerCache{T}) where {T}
    cache.dtacc = one(T)
    cache.erracc = one(T)
    cache.qold = one(T)
    return nothing
end

function sync_controllers!(cache1::PredictiveControllerCache, cache2::PredictiveControllerCache)
    cache1.dtacc = cache2.dtacc
    cache1.erracc = cache2.erracc
    cache1.qold = cache2.qold
    return nothing
end

function setup_controller_cache(alg, cache, controller::PredictiveController, ::Type{E}) where {E}
    QT = _resolved_QT(controller.basic)
    basic = resolve_basic(controller.basic, alg, QT)
    resolved = PredictiveController(basic)
    T = QT
    return PredictiveControllerCache{T, E}(
        resolved, one(T), one(T), one(T), oneunit(E),
    )
end

@inline function stepsize_controller!(integrator, cache::PredictiveControllerCache, alg)
    (; qmin, qmax, gamma) = cache.controller.basic
    qmax = get_current_qmax(integrator, qmax)
    EEst = DiffEqBase.value(get_EEst(integrator))
    if iszero(EEst)
        q = inv(qmax)
    else
        if fac_default_gamma(alg)
            fac = gamma
        else
            if isfirk(alg)
                (; iter) = integrator.cache
                (; maxiters) = alg
            else
                (; iter, maxiters) = integrator.cache.nlsolver
            end
            fac = min(gamma, (1 + 2 * maxiters) * gamma / (iter + 2 * maxiters))
        end
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = fastpower(EEst, expo) / fac
        @fastmath q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        cache.qold = q
    end
    return q
end

function step_accept_controller!(integrator, cache::PredictiveControllerCache, alg, q)
    (; dtacc, erracc, controller) = cache
    (; qmin, qmax, gamma, qsteady_min, qsteady_max) = controller.basic
    qmax = get_current_qmax(integrator, qmax)

    EEst = DiffEqBase.value(get_EEst(integrator))

    if integrator.success_iter > 0
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qgus = (dtacc / integrator.dt) *
            fastpower((EEst^2) / erracc, expo)
        qgus = max(inv(qmax), min(inv(qmin), qgus / gamma))
        qacc = max(q, qgus)
    else
        qacc = q
    end
    if qsteady_min <= qacc <= qsteady_max
        qacc = one(qacc)
    end
    cache.dtacc = DiffEqBase.value(integrator.dt)
    cache.erracc = max(1.0e-2, EEst)

    return integrator.dt / qacc
end

function step_reject_controller!(integrator, cache::PredictiveControllerCache, alg)
    (; dt, success_iter) = integrator
    (; qold) = cache
    return integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
end

# For a composite algorithm the default strategy is to switch forth and back between the controllers of the individual algorithms
struct CompositeController{T} <: AbstractController
    controllers::T
end

mutable struct CompositeControllerCache{T, E} <: AbstractControllerCache
    caches::T
    EEst::E
end

function setup_controller_cache(alg::CompositeAlgorithm, caches::CompositeCache, cc::CompositeController, ::Type{E}) where {E}
    sub = map((alg, cache, controller) -> setup_controller_cache(alg, cache, controller, E), alg.algs, caches.caches, cc.controllers)
    return CompositeControllerCache{typeof(sub), E}(sub, oneunit(E))
end

@inline function accept_step_controller(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return accept_step_controller(integrator, @inbounds(cache.caches[current_idx]), @inbounds(alg.algs[current_idx]))
end
@inline function accept_step_controller(integrator, cache::Union{CompositeCache, CompositeControllerCache}, alg)
    current_idx = integrator.cache.current
    return accept_step_controller(integrator, @inbounds(cache.caches[current_idx]), alg)
end

@inline function stepsize_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return stepsize_controller!(integrator, @inbounds(cache.caches[current_idx]), @inbounds(alg.algs[current_idx]))
end
@inline function stepsize_controller!(integrator, cache::Union{CompositeCache, CompositeControllerCache}, alg)
    current_idx = integrator.cache.current
    return stepsize_controller!(integrator, @inbounds(cache.caches[current_idx]), alg)
end

@inline function step_accept_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm, q)
    current_idx = integrator.cache.current
    return step_accept_controller!(integrator, @inbounds(cache.caches[current_idx]), @inbounds(alg.algs[current_idx]), q)
end
@inline function step_accept_controller!(integrator, cache::Union{CompositeCache, CompositeControllerCache}, alg, q)
    current_idx = integrator.cache.current
    return step_accept_controller!(integrator, @inbounds(cache.caches[current_idx]), alg, q)
end

@inline function step_reject_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return step_reject_controller!(integrator, @inbounds(cache.caches[current_idx]), @inbounds(alg.algs[current_idx]))
end
@inline function step_reject_controller!(integrator, cache::Union{CompositeCache, CompositeControllerCache}, alg)
    current_idx = integrator.cache.current
    return step_reject_controller!(integrator, @inbounds(cache.caches[current_idx]), alg)
end

@inline function post_newton_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return post_newton_controller!(integrator, @inbounds(cache.caches[current_idx]), @inbounds(alg.algs[current_idx]))
end
@inline function post_newton_controller!(integrator, cache::Union{CompositeCache, CompositeControllerCache}, alg)
    current_idx = integrator.cache.current
    return post_newton_controller!(integrator, @inbounds(cache.caches[current_idx]), alg)
end

for accessor in (
        :get_qmin, :get_qmax, :get_qmax_first_step,
        :get_gamma, :get_qsteady_min, :get_qsteady_max,
        :get_failfactor,
    )
    @eval @inline function $accessor(integrator, cache::CompositeControllerCache)
        current_idx = integrator.cache.current
        return $accessor(integrator, @inbounds(cache.caches[current_idx]))
    end
end

function setup_controller_cache(alg::CompositeAlgorithm, caches::DefaultCache, controller::CompositeController, ::Type{E}) where {E}
    sub = (
        setup_controller_cache(alg.algs[1], caches, controller.controllers[1], E),
        setup_controller_cache(alg.algs[2], caches, controller.controllers[2], E),
        setup_controller_cache(alg.algs[3], caches, controller.controllers[3], E),
        setup_controller_cache(alg.algs[4], caches, controller.controllers[4], E),
        setup_controller_cache(alg.algs[5], caches, controller.controllers[5], E),
        setup_controller_cache(alg.algs[6], caches, controller.controllers[6], E),
    )
    return CompositeControllerCache{typeof(sub), E}(sub, oneunit(E))
end

# Default alg
function stepsize_controller!(integrator, cache::DefaultCache, alg)
    return if cache.current == 1
        stepsize_controller!(integrator, @inbounds(cache.cache1), alg)
    elseif cache.current == 2
        stepsize_controller!(integrator, @inbounds(cache.cache2), alg)
    elseif cache.current == 3
        stepsize_controller!(integrator, @inbounds(cache.cache3), alg)
    elseif cache.current == 4
        stepsize_controller!(integrator, @inbounds(cache.cache4), alg)
    elseif cache.current == 5
        stepsize_controller!(integrator, @inbounds(cache.cache5), alg)
    elseif cache.current == 6
        stepsize_controller!(integrator, @inbounds(cache.cache6), alg)
    end
end

function step_accept_controller!(integrator, cache::DefaultCache, alg, q)
    return if cache.current == 1
        step_accept_controller!(integrator, @inbounds(cache.cache1), alg, q)
    elseif cache.current == 2
        step_accept_controller!(integrator, @inbounds(cache.cache2), alg, q)
    elseif cache.current == 3
        step_accept_controller!(integrator, @inbounds(cache.cache3), alg, q)
    elseif cache.current == 4
        step_accept_controller!(integrator, @inbounds(cache.cache4), alg, q)
    elseif cache.current == 5
        step_accept_controller!(integrator, @inbounds(cache.cache5), alg, q)
    elseif cache.current == 6
        step_accept_controller!(integrator, @inbounds(cache.cache6), alg, q)
    end
end

function step_reject_controller!(integrator, cache::DefaultCache, alg)
    return if cache.current == 1
        step_reject_controller!(integrator, @inbounds(cache.cache1), alg)
    elseif cache.current == 2
        step_reject_controller!(integrator, @inbounds(cache.cache2), alg)
    elseif cache.current == 3
        step_reject_controller!(integrator, @inbounds(cache.cache3), alg)
    elseif cache.current == 4
        step_reject_controller!(integrator, @inbounds(cache.cache4), alg)
    elseif cache.current == 5
        step_reject_controller!(integrator, @inbounds(cache.cache5), alg)
    elseif cache.current == 6
        step_reject_controller!(integrator, @inbounds(cache.cache6), alg)
    end
end

# This is a workaround to make the BDF methods work with composite algorithms
function post_newton_controller!(integrator, cache::DefaultCache, alg)
    if cache.current == 1
        post_newton_controller!(integrator, @inbounds(cache.cache1), alg)
    elseif cache.current == 2
        post_newton_controller!(integrator, @inbounds(cache.cache2), alg)
    elseif cache.current == 3
        post_newton_controller!(integrator, @inbounds(cache.cache3), alg)
    elseif cache.current == 4
        post_newton_controller!(integrator, @inbounds(cache.cache4), alg)
    elseif cache.current == 5
        post_newton_controller!(integrator, @inbounds(cache.cache5), alg)
    elseif cache.current == 6
        post_newton_controller!(integrator, @inbounds(cache.cache6), alg)
    end
    return nothing
end
