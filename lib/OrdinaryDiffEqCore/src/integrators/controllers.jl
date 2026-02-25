abstract type AbstractController end
abstract type AbstractLegacyController <: AbstractController end

abstract type AbstractControllerCache end

"""
    setup_controller_cache(alg, atmp, controller::AbstractController)::AbstractControllerCache

This function takes a controller together with the time stepping algorithm to
construct and initialize the respective cache for the controller.
"""
setup_controller_cache

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

!!! warning
    The parameter `q` will be removed the next release and will be passed through the
    controller cache if needed.
"""
step_accept_controller!

"""
    step_reject_controller!(integrator, alg)
    step_reject_controller!(integrator, controller_cache::AbstractControllerCache, alg)

This function gets called in case of a rejected time step right after [`stepsize_controller!`](@ref).
It directly sets the time step length (i.e. `integrator.dt`).
"""
step_reject_controller!


# The legacy controllers do not have this concept.
setup_controller_cache(alg, atmp, controller::AbstractLegacyController) = controller

# checks whether the controller should accept a step based on the error estimate
@inline function accept_step_controller(integrator, controller_or_cache::Union{<:AbstractLegacyController, <:AbstractControllerCache})
    return integrator.EEst <= 1
end

@inline function stepsize_controller!(integrator, alg)
    # TODO replace this when done - right now this holds the controller cache!
    return stepsize_controller!(integrator, integrator.opts.controller, alg)
    # stepsize_controller!(integrator, integrator.controller_cache, alg)
    # return nothing
end

@inline function step_accept_controller!(integrator, alg, q)
    # TODO replace this when done - right now this holds the controller cache!
    return step_accept_controller!(integrator, integrator.opts.controller, alg, q)
    # step_accept_controller!(integrator, integrator.controller_cache, alg, q)
    # return nothing
end

@inline function step_reject_controller!(integrator, alg)
    # TODO replace this when done - right now this holds the controller cache!
    step_reject_controller!(integrator, integrator.opts.controller, alg)
    # step_reject_controller!(integrator, integrator.controller_cache, alg)
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
        ctrl = integrator.opts.controller
        if hasfield(typeof(ctrl), :controller) &&
                hasfield(typeof(ctrl.controller), :qmax_first_step)
            return ctrl.controller.qmax_first_step
        end
        return typeof(qmax)(10000)
    end
    return qmax
end

reset_alg_dependent_opts!(controller::AbstractController, alg1, alg2) = nothing
reset_alg_dependent_opts!(controller::AbstractControllerCache, alg1, alg2) = nothing

SciMLBase.reinit!(integrator::ODEIntegrator, controller::AbstractController) = nothing
SciMLBase.reinit!(integrator::ODEIntegrator, controller::AbstractControllerCache) = nothing

function post_newton_controller!(integrator, alg)
    return post_newton_controller!(integrator, integrator.opts.controller, alg)
end
function post_newton_controller!(integrator, controller, alg)
    integrator.dt = integrator.dt / integrator.opts.failfactor
    return nothing
end

# This is a helper for Nordsieck and BDF methods, which come with an integrated controller.
struct DummyController <: AbstractController
end

setup_controller_cache(alg, atmp, controller::DummyController) = controller

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
inverse of the error estimate `integrator.EEst` scaled by the tolerance
(Hairer, Nørsett, Wanner, 2008, Section II.4).
The step size factor is multiplied by the safety factor `gamma` and clipped to
the interval `[qmin, qmax]`.
A step will be accepted whenever the estimated error `integrator.EEst` is
less than or equal to unity. Otherwise, the step is rejected and re-tried with
the predicted step size.

## References

  - Hairer, Nørsett, Wanner (2008)
    Solving Ordinary Differential Equations I Nonstiff Problems
    [DOI: 10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1)
"""
struct IController <: AbstractLegacyController
end

@inline function stepsize_controller!(integrator, controller::IController, alg)
    (; qmin, qmax, gamma) = integrator.opts
    qmax = get_current_qmax(integrator, qmax)
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = inv(qmax)
    else
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = fastpower(EEst, expo) / gamma
        @fastmath q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        # TODO: Shouldn't this be in `step_accept_controller!` as for the PI controller?
        integrator.qold = DiffEqBase.value(integrator.dt) / q
    end
    return q
end

function step_accept_controller!(integrator, controller::IController, alg, q)
    (; qsteady_min, qsteady_max) = integrator.opts

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, controller::IController, alg)
    (; qold) = integrator
    return integrator.dt = qold
end

struct NewIController{T} <: AbstractController
    qmin::T
    qmax::T
    qmax_first_step::T
    gamma::T
    qsteady_min::T
    qsteady_max::T
end

function NewIController(alg; kwargs...)
    return NewIController(Float64, alg; kwargs...)
end

function NewIController(QT, alg; qmin = nothing, qmax = nothing, qmax_first_step = nothing, gamma = nothing, qsteady_min = nothing, qsteady_max = nothing)
    return NewIController{QT}(
        qmin === nothing ? qmin_default(alg) : qmin,
        qmax === nothing ? qmax_default(alg) : qmax,
        qmax_first_step === nothing ? QT(10000) : QT(qmax_first_step),
        gamma === nothing ? gamma_default(alg) : gamma,
        qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min,
        qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max,
    )
end

mutable struct IControllerCache{C, T, UT} <: AbstractControllerCache
    controller::C
    q::T
    dtreject::T
    # I believe this should go here or in the algorithm cache, but not in the integrator itself.
    # EEst::T
    atmp::UT
end

function setup_controller_cache(alg, atmp, controller::NewIController{T}) where {T}
    return IControllerCache(
        controller,
        one(T),
        T(1 // 10^4), # TODO which value?
        atmp,
    )
end

@inline function stepsize_controller!(integrator, cache::IControllerCache, alg)
    (; qmin, qmax, gamma) = cache.controller
    qmax = get_current_qmax(integrator, qmax)
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        cache.q = inv(qmax)
    else
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = fastpower(EEst, expo) / gamma
        @fastmath cache.q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        # TODO: Shouldn't this be in `step_accept_controller!` as for the PI controller?
        cache.dtreject = integrator.qold = DiffEqBase.value(integrator.dt) / cache.q
    end
    return cache.q
end

# TODO change signature to remove the q input
function step_accept_controller!(integrator, cache::IControllerCache, alg, q)
    (; qsteady_min, qsteady_max) = cache.controller
    @assert q ≈ cache.q "Controller cache went out of sync with time stepping logic."

    if qsteady_min <= q <= qsteady_max
        cache.q = q = one(q)
    end
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, cache::IControllerCache, alg)
    @assert cache.dtreject ≈ integrator.qold "Controller cache went out of sync with time stepping logic."
    return integrator.dt = cache.dtreject # TODO this does not look right.
end


SciMLBase.reinit!(integrator::ODEIntegrator, cache::IControllerCache{T}) where {T} = cache.q = one(T)

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
A step will be accepted whenever the estimated error `integrator.EEst` is
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
mutable struct PIController{QT} <: AbstractLegacyController
    beta1::QT
    beta2::QT
end

@inline function stepsize_controller!(integrator, controller::PIController, alg)
    (; qold) = integrator
    (; qmin, qmax, gamma) = integrator.opts
    qmax = get_current_qmax(integrator, qmax)
    (; beta1, beta2) = controller
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = inv(qmax)
    else
        q11 = fastpower(EEst, convert(typeof(EEst), beta1))
        q = q11 / fastpower(qold, convert(typeof(EEst), beta2))
        integrator.q11 = q11
        @fastmath q = max(inv(qmax), min(inv(qmin), q / gamma))
    end
    return q
end

function step_accept_controller!(integrator, controller::PIController, alg, q)
    (; qsteady_min, qsteady_max, qoldinit) = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    integrator.qold = max(EEst, qoldinit)
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, controller::PIController, alg)
    (; q11) = integrator
    (; qmin, gamma) = integrator.opts
    return integrator.dt /= min(inv(qmin), q11 / gamma)
end

function reset_alg_dependent_opts!(controller::PIController{QT}, alg1, alg2) where {QT}
    # These guards here prevent potentially that user-defined betas get overridden.
    if controller.beta2 == QT(beta2_default(alg1))
        controller.beta2 = QT(beta2_default(alg2))
    end
    if controller.beta1 == QT(beta1_default(alg1, controller.beta2))
        controller.beta1 = QT(beta1_default(alg2, controller.beta2))
    end
    return nothing
end

struct NewPIController{T} <: AbstractController
    beta1::T
    beta2::T
    qmin::T
    qmax::T
    qmax_first_step::T
    gamma::T
    qsteady_min::T
    qsteady_max::T
    qoldinit::T
end

function NewPIController(alg; kwargs...)
    return NewPIController(Float64, alg; kwargs...)
end

function NewPIController(QT, alg; beta1 = nothing, beta2 = nothing, qmin = nothing, qmax = nothing, qmax_first_step = nothing, gamma = nothing, qsteady_min = nothing, qsteady_max = nothing, qoldinit = nothing)
    beta2 = beta2 === nothing ? beta2_default(alg) : beta2
    beta1 = beta1 === nothing ? beta1_default(alg, beta2) : beta1
    qoldinit = qoldinit === nothing ? 1 // 10^4 : qoldinit
    return NewPIController{QT}(
        beta1,
        beta2,
        qmin === nothing ? qmin_default(alg) : qmin,
        qmax === nothing ? qmax_default(alg) : qmax,
        qmax_first_step === nothing ? QT(10000) : QT(qmax_first_step),
        gamma === nothing ? gamma_default(alg) : gamma,
        qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min,
        qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max,
        qoldinit,
    )
end

mutable struct PIControllerCache{T, UT} <: AbstractControllerCache
    controller::NewPIController{T}
    # Propsoed scaling factor for the time step length
    q::T
    # Cached εₙ₊₁^β₁
    q11::T
    # Previous EEst
    errold::T
    atmp::UT
end

function setup_controller_cache(alg, atmp, controller::NewPIController{T}) where {T}
    return PIControllerCache(
        controller,
        one(T),
        one(T),
        T(controller.qoldinit),
        atmp,
    )
end

@inline function stepsize_controller!(integrator, cache::PIControllerCache, alg)
    (; errold, controller) = cache
    (; qmin, qmax, gamma) = controller
    qmax = get_current_qmax(integrator, qmax)
    (; beta1, beta2) = controller
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = inv(qmax)
    else
        # Legacy code
        q11 = fastpower(EEst, beta1)
        q = q11 / fastpower(errold, beta2)
        cache.q11 = q11
        integrator.q11 = q11 # TODO remove
        @fastmath q = clamp(q / gamma, inv(qmax), inv(qmin))
    end
    cache.q = q
    return q
end

function step_accept_controller!(integrator, cache::PIControllerCache, alg, q)
    @assert q ≈ cache.q "Controller cache went out of sync with time stepping logic (q=$q | cache.q=$(cache.q))."
    (; controller) = cache
    (; qsteady_min, qsteady_max, qoldinit) = controller
    EEst = DiffEqBase.value(integrator.EEst)

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    cache.errold = max(EEst, qoldinit)
    integrator.qold = cache.errold # TODO remove
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, cache::PIControllerCache, alg)
    (; controller, q11) = cache
    (; qmin, gamma) = controller
    return integrator.dt /= min(inv(qmin), q11 / gamma)
end

function SciMLBase.reinit!(integrator::ODEIntegrator, cache::PIControllerCache{T}) where {T}
    cache.q = one(T)
    cache.q11 = one(T)
    return cache.errold = T(cache.controller.qoldinit)
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
    [DOI: 10.1016/j.cam.2005.03.008](https://doi.org/10.1016/j.cam.2005.03.008) # controller coefficients
  - Ranocha, Dalcin, Parsani, Ketcheson (2021) # history of the error estimates
    Optimized Runge-Kutta Methods with Automatic Step Size Control for   # accept a step if the predicted change of the step size
    Compressible Computational Fluid Dynamics    # is bigger than this parameter
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)    # limiter of the dt factor (before clipping)
"""
struct PIDController{QT, Limiter} <: AbstractLegacyController
    beta::MVector{3, QT} # controller coefficients
    err::MVector{3, QT} # history of the error estimates
    accept_safety::QT   # accept a step if the predicted change of the step size
    # is bigger than this parameter
    limiter::Limiter    # limiter of the dt factor (before clipping)
end

function PIDController(
        beta1, beta2, beta3 = zero(beta1);
        limiter = default_dt_factor_limiter,
        accept_safety = 0.81
    )
    beta = MVector(map(float, promote(beta1, beta2, beta3))...)
    QT = eltype(beta)
    err = MVector{3, QT}(true, true, true)
    return PIDController(beta, err, convert(QT, accept_safety), limiter)
end

function Base.show(io::IO, controller::PIDController)
    return print(
        io, "PIDController(beta=", controller.beta,
        ", accept_safety=", controller.accept_safety,
        ", limiter=", controller.limiter,
        ")"
    )
end

@inline default_dt_factor_limiter(x) = one(x) + atan(x - one(x))

@inline function stepsize_controller!(integrator, controller::PIDController, alg)
    beta1, beta2, beta3 = controller.beta

    EEst = DiffEqBase.value(integrator.EEst)

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

    controller.err[1] = inv(EEst)
    err1, err2, err3 = controller.err

    k = min(alg_order(alg), alg_adaptive_order(alg)) + 1
    dt_factor = err1^(beta1 / k) * err2^(beta2 / k) * err3^(beta3 / k)
    if isnan(dt_factor)
        @SciMLMessage(
            lazy"unlimited dt_factor:
            err1 = $err1
            err2 = $err2,
            err3 = $err3,
            beta1 = $beta1,
            beta2 = $beta2,
            beta3 = $beta3,
            k = $k",
            integrator.opts.verbose, :unlimited_dt
        )
    end
    dt_factor = controller.limiter(dt_factor)

    # Note: No additional limiting of the form
    #   dt_factor = max(qmin, min(qmax, dt_factor))
    # is necessary since the `limiter` should take care of that. The default limiter
    # ensures
    #   0.21 ≈ limiter(0) <= dt_factor <= limiter(Inf) ≈ 2.57
    # See Söderlind, Wang (2006), Section 6.
    integrator.qold = dt_factor
    return dt_factor
end

@inline function accept_step_controller(integrator, controller::PIDController)
    return integrator.qold >= controller.accept_safety
end

function step_accept_controller!(integrator, controller::PIDController, alg, dt_factor)
    (; qsteady_min, qsteady_max) = integrator.opts

    if qsteady_min <= inv(dt_factor) <= qsteady_max
        dt_factor = one(dt_factor)
    end
    @inbounds begin
        controller.err[3] = controller.err[2]
        controller.err[2] = controller.err[1]
    end
    return integrator.dt * dt_factor # new dt
end

function step_reject_controller!(integrator, controller::PIDController, alg)
    return integrator.dt *= integrator.qold
end


struct NewPIDController{T, Limiter} <: AbstractController
    beta::SVector{3, T} # controller coefficients
    accept_safety::T   # accept a step if the predicted change of the step size
    # is bigger than this parameter
    limiter::Limiter    # limiter of the dt factor (before clipping)
    qsteady_min::T
    qsteady_max::T
end

function NewPIDController(alg; kwargs...)
    return NewPIDController(Float64, alg; kwargs...)
end

function NewPIDController(QT, alg; beta = nothing, accept_safety = 0.81, limiter = default_dt_factor_limiter, qsteady_min = nothing, qsteady_max = nothing)
    if beta === nothing
        beta2 = QT(beta2_default(alg))
        beta1 = QT(beta1_default(alg, beta2))
        beta3 = QT(zero(beta1))
    else
        beta1, beta2, beta3 = beta
    end
    beta = SVector(map(float, promote(beta1, beta2, beta3))...)
    return NewPIDController{QT, typeof(limiter)}(
        beta,
        QT(accept_safety),
        limiter,
        QT(qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min),
        QT(qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max),
    )
end

function Base.show(io::IO, controller::NewPIDController)
    return print(
        io, "NewPIDController(beta=", controller.beta,
        ", accept_safety=", controller.accept_safety,
        ", limiter=", controller.limiter,
        ")"
    )
end

mutable struct PIDControllerCache{T, Limiter, UT} <: AbstractControllerCache
    controller::NewPIDController{T, Limiter}
    err::MVector{3, T} # history of the error estimates
    dt_factor::T
    atmp::UT
end

function SciMLBase.reinit!(integrator::ODEIntegrator, cache::PIDControllerCache{T}) where {T}
    cache.err = MVector{3, T}(true, true, true)
    return cache.dt_factor = T(1 // 10^4)
end

function setup_controller_cache(alg, atmp, controller::NewPIDController{QT}) where {QT}
    err = MVector{3, QT}(true, true, true)
    return PIDControllerCache(
        controller,
        err,
        QT(1 // 10^4),
        atmp,
    )
end

@inline function stepsize_controller!(integrator, cache::PIDControllerCache, alg)
    (; controller) = cache
    beta1, beta2, beta3 = controller.beta

    EEst = DiffEqBase.value(integrator.EEst)

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
    dt_factor = controller.limiter(dt_factor)

    # Note: No additional limiting of the form
    #   dt_factor = max(qmin, min(qmax, dt_factor))
    # is necessary since the `limiter` should take care of that. The default limiter
    # ensures
    #   0.21 ≈ limiter(0) <= dt_factor <= limiter(Inf) ≈ 2.57
    # See Söderlind, Wang (2006), Section 6.
    cache.dt_factor = dt_factor
    return dt_factor
end

@inline function accept_step_controller(integrator, cache::PIDControllerCache)
    return cache.dt_factor >= cache.controller.accept_safety
end

function step_accept_controller!(integrator, cache::PIDControllerCache, alg, dt_factor)
    @assert dt_factor ≈ cache.dt_factor "Controller cache went out of sync with time stepping logic."
    (; controller) = cache
    (; qsteady_min, qsteady_max) = controller

    if qsteady_min <= inv(dt_factor) <= qsteady_max
        dt_factor = one(dt_factor)
        # cache.q = ...?
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

# Gustafsson predictive step size controller
"""
    PredictiveController()

The Gustafsson acceleration algorithm accelerates changes so that way algorithms
can more swiftly change to handle quick transients. This algorithm is thus
well-suited for stiff solvers where this can be expected, and is the default
for algorithms like the (E)SDIRK methods.

```julia
gamma = integrator.opts.gamma
niters = integrator.cache.newton_iters
fac = min(gamma,
    (1 + 2 * integrator.alg.max_newton_iter) * gamma /
    (niters + 2 * integrator.alg.max_newton_iter))
expo = 1 / (alg_order(integrator.alg) + 1)
qtmp = (integrator.EEst^expo) / fac
@fastmath q = max(inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), qtmp))
if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
    q = one(q)
end
integrator.qold = q
q
```

In this case, `niters` is the number of Newton iterations which was required
in the most recent step of the algorithm. Note that these values are used
differently depending on acceptance and rejectance. When the step is accepted,
the following logic is applied:

```julia
if integrator.success_iter > 0
    expo = 1 / (alg_adaptive_order(integrator.alg) + 1)
    qgus = (integrator.dtacc / integrator.dt) * (((integrator.EEst^2) / integrator.erracc)^expo)
    qgus = max(inv(integrator.opts.qmax),
        min(inv(integrator.opts.qmin), qgus / integrator.opts.gamma))
    qacc = max(q, qgus)
else
    qacc = q
end
integrator.dtacc = integrator.dt
integrator.erracc = max(1e-2, integrator.EEst)
integrator.dt / qacc
```

When it rejects, it's the same as the [`IController`](@ref):

```julia
if integrator.success_iter == 0
    integrator.dt *= 0.1
else
    integrator.dt = integrator.dt / integrator.qold
end
```
"""
struct PredictiveController <: AbstractLegacyController
end

@inline function stepsize_controller!(integrator, controller::PredictiveController, alg)
    (; qmin, qmax, gamma) = integrator.opts
    qmax = get_current_qmax(integrator, qmax)
    EEst = DiffEqBase.value(integrator.EEst)
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
        integrator.qold = q
    end
    return q
end

function step_accept_controller!(integrator, controller::PredictiveController, alg, q)
    (; qmin, qmax, gamma, qsteady_min, qsteady_max) = integrator.opts
    qmax = get_current_qmax(integrator, qmax)

    EEst = DiffEqBase.value(integrator.EEst)

    if integrator.success_iter > 0
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qgus = (integrator.dtacc / integrator.dt) *
            fastpower((EEst^2) / integrator.erracc, expo)
        qgus = max(inv(qmax), min(inv(qmin), qgus / gamma))
        qacc = max(q, qgus)
    else
        qacc = q
    end
    if qsteady_min <= qacc <= qsteady_max
        qacc = one(qacc)
    end
    integrator.dtacc = integrator.dt
    integrator.erracc = max(1.0e-2, EEst)

    return integrator.dt / qacc
end

function step_reject_controller!(integrator, controller::PredictiveController, alg)
    (; dt, success_iter, qold) = integrator
    return integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
end


struct NewPredictiveController{T} <: AbstractController
    qmin::T
    qmax::T
    qmax_first_step::T
    gamma::T
    qsteady_min::T
    qsteady_max::T
end

function NewPredictiveController(alg; kwargs...)
    return NewPredictiveController(Float64, alg; kwargs...)
end

function NewPredictiveController(QT, alg; qmin = nothing, qmax = nothing, qmax_first_step = nothing, gamma = nothing, qsteady_min = nothing, qsteady_max = nothing)
    return NewPredictiveController{QT}(
        qmin === nothing ? qmin_default(alg) : qmin,
        qmax === nothing ? qmax_default(alg) : qmax,
        qmax_first_step === nothing ? QT(10000) : QT(qmax_first_step),
        gamma === nothing ? gamma_default(alg) : gamma,
        qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min,
        qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max,
    )
end

mutable struct PredictiveControllerCache{T, UT} <: AbstractControllerCache
    controller::NewPredictiveController{T}
    dtacc::T
    erracc::T
    qold::T
    q::T
    atmp::UT
end

function SciMLBase.reinit!(integrator::ODEIntegrator, cache::PredictiveControllerCache{T}) where {T}
    cache.dtacc = one(T)
    cache.erracc = one(T)
    cache.qold = one(T)
    return cache.q = one(T)
end

function setup_controller_cache(alg, atmp::UT, controller::NewPredictiveController{T}) where {T, UT}
    return PredictiveControllerCache{T, UT}(
        controller,
        one(T),
        one(T),
        one(T),
        one(T),
        atmp,
    )
end

@inline function stepsize_controller!(integrator, cache::PredictiveControllerCache, alg)
    (; qmin, qmax, gamma) = cache.controller
    qmax = get_current_qmax(integrator, qmax)
    EEst = DiffEqBase.value(integrator.EEst)
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
    cache.q = q
    return q
end

function step_accept_controller!(integrator, cache::PredictiveControllerCache, alg, q)
    @assert q ≈ cache.q "Controller cache went out of sync with time stepping logic."
    (; dtacc, erracc, controller) = cache
    (; qmin, qmax, gamma, qsteady_min, qsteady_max) = controller
    qmax = get_current_qmax(integrator, qmax)

    EEst = DiffEqBase.value(integrator.EEst)

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

struct CompositeControllerCache{T, UT} <: AbstractControllerCache
    caches::T
    atmp::UT # This is just here for easy access
end

function setup_controller_cache(alg::CompositeAlgorithm, atmp, cc::CompositeController)
    return CompositeControllerCache(
        map((alg, controller) -> setup_controller_cache(alg, atmp, controller), alg.algs, cc.controllers),
        atmp,
    )
end

@inline function accept_step_controller(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return accept_step_controller(integrator, cache.caches[current_idx], alg.algs[current_idx])
end

@inline function stepsize_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return stepsize_controller!(integrator, cache.caches[current_idx], alg.algs[current_idx])
end

@inline function step_accept_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm, q)
    current_idx = integrator.cache.current
    return step_accept_controller!(integrator, cache.caches[current_idx], alg.algs[current_idx], q)
end

@inline function step_reject_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return step_reject_controller!(integrator, cache.caches[current_idx], alg.algs[current_idx])
end

@inline function post_newton_controller!(integrator, cache::CompositeControllerCache, alg::CompositeAlgorithm)
    current_idx = integrator.cache.current
    return post_newton_controller!(integrator, cache.caches[current_idx], alg.algs[current_idx])
end

# We need this for now as a workaround to make composite controllers work when they have a dummy controller attached
@inline function accept_step_controller(integrator, controller::DummyController)
    return integrator.EEst <= 1
end
@inline stepsize_controller!(integrator, controller::DummyController, alg) = default_stepsize_controller!(integrator, integrator.cache, alg)
@inline step_accept_controller!(integrator, controller::DummyController, alg, q) = default_step_accept_controller!(integrator, integrator.cache, alg, q)
@inline step_reject_controller!(integrator, controller::DummyController, alg) = default_step_reject_controller!(integrator, integrator.cache, alg)
@inline post_newton_controller!(integrator, controller::DummyController, alg) = default_post_newton_controller!(integrator, integrator.cache, alg)

# TODO remove this for OrdinaryDiffEq v7 . Right now the integrator is expected to carry a controller. Therefore algorithms coming with a custom controller default to a DummyController too.
# Instead of scattering this function across all subpackages we add the default dispatch here.
function default_post_newton_controller!(integrator, cache, alg)
    integrator.dt = integrator.dt / integrator.opts.failfactor
    return nothing
end

# Default alg with dummy controller
function default_stepsize_controller!(integrator, cache::DefaultCache, alg)
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

function default_step_accept_controller!(integrator, cache::DefaultCache, alg, q)
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

function default_step_reject_controller!(integrator, cache::DefaultCache, alg)
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

# Composite alg with dummy controller
function default_stepsize_controller!(integrator, cache::CompositeCache, alg)
    return stepsize_controller!(integrator, @inbounds(cache.caches[cache.current]), alg)
end

function default_step_accept_controller!(integrator, cache::CompositeCache, alg, q)
    return step_accept_controller!(integrator, @inbounds(cache.caches[cache.current]), alg, q)
end

function default_step_reject_controller!(integrator, cache::CompositeCache, alg)
    return step_reject_controller!(integrator, @inbounds(cache.caches[cache.current]), alg)
end

# This is a workaround to make the BDF methods work with composite algorithms
function default_post_newton_controller!(integrator, cache::DefaultCache, alg)
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
function default_post_newton_controller!(integrator, cache::CompositeCache, alg)
    post_newton_controller!(integrator, @inbounds(cache.caches[cache.current]), alg)
    return nothing
end
