abstract type AbstractController end
abstract type AbstractLegacyController <: AbstractController end

abstract type AbstractControllerCache end

"""
    setup_controller_cache(alg_cache, controller::AbstractController)::AbstractControllerCache

This function takes a controller together with the cache of the time stepping
algorithm to construct and initialize the respective cache for the controller.
"""
setup_controller_cache

"""
    accept_step_controller(integrator, controller_cache::AbstractControllerCache)::Bool

This function decides whether the current time step should be accepted or rejected.
A return value of `false` corresponds to a rejection.
"""
accept_step_controller

"""
    stepsize_controller!(integrator, controller_cache::AbstractControllerCache)

Update the cache to compute the new order of the time marching algorithm and prepare
for an update of the time step. The update can be either due to a rejection or acceptance
of the current time step.
"""
stepsize_controller!

"""
    step_accept_controller!(integrator, controller_cache::AbstractControllerCache, q)

This function gets called in case of an accepted time step right after [`stepsize_controller!`](@ref).
It returns the proposed new time step length. Please note that the time step length might not be
applied as is and subject to further modification to e.g. match the next time stop.

!!! warning
    The parameter `q` will be removed the next release and will be passed through the
    controller cache if needed.
"""
step_accept_controller!

"""
    step_reject_controller!(integrator, controller_cache::AbstractControllerCache)

This function gets called in case of a rejected time step right after [`stepsize_controller!`](@ref).
It directly sets the time step length (i.e. `integrator.dt`).
"""
step_reject_controller!


# The legacy controllers do not have this concept.
setup_controller_cache(alg_cache, controller::AbstractLegacyController) = controller

# checks whether the controller should accept a step based on the error estimate
@inline function accept_step_controller(integrator, controller_or_cache::Union{<:AbstractLegacyController, <:AbstractControllerCache})
    return integrator.EEst <= 1
end

@inline function stepsize_controller!(integrator, alg)
    stepsize_controller!(integrator, integrator.controller_cache, alg)
end

@inline function step_accept_controller!(integrator, alg, q)
    step_accept_controller!(integrator, integrator.controller_cache, alg, q)
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

reset_alg_dependent_opts!(controller::AbstractController, alg1, alg2) = nothing

SciMLBase.reinit!(integrator::ODEIntegrator, controller::AbstractController) = nothing

function post_newton_controller!(integrator, alg)
    integrator.dt = integrator.dt / integrator.opts.failfactor
    nothing
end

# This is a helper for Nordsieck and BDF methods, which come with an integrated controller.
struct DummyController <: AbstractController
end

setup_controller_cache(alg_cache, controller::DummyController) = controller

@inline function accept_step_controller(integrator, controller::DummyController)
    return integrator.EEst <= 1
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
    @unpack qmin, qmax, gamma = integrator.opts
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
    q
end

function step_accept_controller!(integrator, controller::IController, alg, q)
    @unpack qsteady_min, qsteady_max = integrator.opts

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    integrator.dt / q # new dt
end

function step_reject_controller!(integrator, controller::IController, alg)
    @unpack qold = integrator
    integrator.dt = qold
end

struct NewIController{T} <: AbstractController
    qmin::T
    qmax::T
    gamma::T
    qsteady_min::T
    qsteady_max::T
end

function NewIController(alg; kwargs...)
    NewIController(Float64, alg; kwargs...)
end

function NewIController(QT, alg; qmin = nothing, qmax = nothing, gamma = nothing, qsteady_min = nothing, qsteady_max = nothing)
    NewIController{QT}(
        qmin === nothing ? qmin_default(alg) : qmin,
        qmax === nothing ? qmax_default(alg) : qmax,
        gamma === nothing ? gamma_default(alg) : gamma,
        qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min,
        qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max,
    )
end

mutable struct IControllerCache{C, T} <: AbstractControllerCache
    controller::C
    q::T
    dtreject::T
    # I believe this should go here or in the algorithm cache, but not in the integrator itself.
    # EEst::T
end

function setup_controller_cache(alg_cache, controller::NewIController{T}) where T
    IControllerCache(
        controller,
        T(1),
        T(1 // 10^4), # TODO which value?
    )
end

@inline function stepsize_controller!(integrator, cache::IControllerCache, alg)
    @unpack qmin, qmax, gamma = cache.controller
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
    cache.q
end

# TODO change signature to remove the q input
function step_accept_controller!(integrator, cache::IControllerCache, alg, q)
    @unpack qsteady_min, qsteady_max = cache.controller
    @assert q ≈ cache.q "Controller cache went out of sync with time stepping logic."

    if qsteady_min <= q <= qsteady_max
        cache.q = q = one(q)
    end
    integrator.dt / q # new dt
end

function step_reject_controller!(integrator, cache::IControllerCache, alg)
    @assert cache.dtreject ≈ integrator.qold "Controller cache went out of sync with time stepping logic."
    integrator.dt = cache.dtreject # TODO this does not look right.
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
    @unpack qold = integrator
    @unpack qmin, qmax, gamma = integrator.opts
    @unpack beta1, beta2 = controller
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = inv(qmax)
    else
        q11 = fastpower(EEst, convert(typeof(EEst), beta1))
        q = q11 / fastpower(qold, convert(typeof(EEst), beta2))
        integrator.q11 = q11
        @fastmath q = max(inv(qmax), min(inv(qmin), q / gamma))
    end
    q
end

function step_accept_controller!(integrator, controller::PIController, alg, q)
    @unpack qsteady_min, qsteady_max, qoldinit = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    integrator.qold = max(EEst, qoldinit)
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, controller::PIController, alg)
    @unpack q11 = integrator
    @unpack qmin, gamma = integrator.opts
    integrator.dt /= min(inv(qmin), q11 / gamma)
end

function reset_alg_dependent_opts!(controller::PIController, alg1, alg2)
    if controller.beta2 == beta2_default(alg1)
        controller.beta2 = beta2_default(alg2)
    end
    if controller.beta1 == beta1_default(alg1, controller.beta2)
        controller.beta1 = beta1_default(alg2, controller.beta2)
    end
end

struct NewPIController{T} <: AbstractController
    beta1::T
    beta2::T
    qmin::T
    qmax::T
    gamma::T
    qsteady_min::T
    qsteady_max::T
    qoldinit::T
end

function NewPIController(alg; kwargs...)
    NewPIController(Float64, alg; kwargs...)
end

function NewPIController(QT, alg; beta1 = nothing, beta2 = nothing, qmin = nothing, qmax = nothing, gamma = nothing, qsteady_min = nothing, qsteady_max = nothing, qoldinit = nothing)
    beta2 = beta2 === nothing ? beta2_default(alg) : beta2
    beta1 = beta1 === nothing ? beta1_default(alg, beta2) : beta1
    qoldinit = qoldinit === nothing ? 1 // 10^4 : qoldinit
    NewPIController{QT}(
        beta1,
        beta2,
        qmin === nothing ? qmin_default(alg) : qmin,
        qmax === nothing ? qmax_default(alg) : qmax,
        gamma === nothing ? gamma_default(alg) : gamma,
        qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min,
        qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max,
        qoldinit,
    )
end

mutable struct PIControllerCache{T} <: AbstractControllerCache
    controller::NewPIController{T}
    # Propsoed scaling factor for the time step length
    q::T
    # Cached εₙ₊₁^β₁
    q11::T
    # Previous EEst
    errold::T
end


function setup_controller_cache(alg_cache, controller::NewPIController{T}) where T
    PIControllerCache(
        controller,
        T(1),
        T(1),
        T(1 // 10^4),
    )
end

@inline function stepsize_controller!(integrator, cache::PIControllerCache, alg)
    @unpack errold, controller = cache
    @unpack qmin, qmax, gamma = controller
    @unpack beta1, beta2 = controller
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = qmax
    else
        # Legacy code
        q11 = fastpower(EEst, beta1)
        q = q11 / fastpower(errold, beta2)
        cache.q11 = q11
        integrator.q11 = q11 # TODO remove
        @fastmath q = clamp(q / gamma, inv(qmax), inv(qmin))
    end
    cache.q = q
    q
end

function step_accept_controller!(integrator, cache::PIControllerCache, alg, q)
    @assert q ≈ cache.q "Controller cache went out of sync with time stepping logic (q=$q | cache.q=$(cache.q))."
    @unpack controller = cache
    @unpack qsteady_min, qsteady_max, qoldinit = controller
    EEst = DiffEqBase.value(integrator.EEst)

    if qsteady_min <= q <= qsteady_max
        q = one(q)
    end
    cache.errold = max(EEst, qoldinit)
    integrator.qold = cache.errold # TODO remove
    return integrator.dt / q # new dt
end

function step_reject_controller!(integrator, cache::PIControllerCache, alg)
    @unpack controller, q11 = cache
    @unpack qmin, gamma = controller
    integrator.dt /= min(inv(qmin), q11 / gamma)
end

# FIXME multi-controller?
# function reset_alg_dependent_opts!(controller::PIControllerCache, alg1, alg2)
#     if controller.beta2 == beta2_default(alg1)
#         controller.beta2 = beta2_default(alg2)
#     end
#     if controller.beta1 == beta1_default(alg1, controller.beta2)
#         controller.beta1 = beta1_default(alg2, controller.beta2)
#     end
# end

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

function PIDController(beta1, beta2, beta3 = zero(beta1);
        limiter = default_dt_factor_limiter,
        accept_safety = 0.81)
    beta = MVector(map(float, promote(beta1, beta2, beta3))...)
    QT = eltype(beta)
    err = MVector{3, QT}(true, true, true)
    return PIDController(beta, err, convert(QT, accept_safety), limiter)
end

function Base.show(io::IO, controller::PIDController)
    print(io, "PIDController(beta=", controller.beta,
        ", accept_safety=", controller.accept_safety,
        ", limiter=", controller.limiter,
        ")")
end

@inline default_dt_factor_limiter(x) = one(x) + atan(x - one(x))

@inline function stepsize_controller!(integrator, controller::PIDController, alg)
    @unpack qmax = integrator.opts
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
        @warn "unlimited dt_factor" dt_factor err1 err2 err3 beta1 beta2 beta3 k
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
    @unpack qsteady_min, qsteady_max = integrator.opts

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
    integrator.dt *= integrator.qold
end



struct NewPIDController{T, Limiter} <: AbstractController
    beta::SVector{3, T} # controller coefficients
    accept_safety::T   # accept a step if the predicted change of the step size
    # is bigger than this parameter
    limiter::Limiter    # limiter of the dt factor (before clipping)
    qmin::T
    qmax::T
    qsteady_min::T
    qsteady_max::T
end

function NewPIDController(alg; kwargs...)
    NewPIDController(Float64, alg; kwargs...)
end

function NewPIDController(QT, alg; beta = nothing, accept_safety = 0.81, limiter = default_dt_factor_limiter, qmin = nothing, qmax = nothing, qsteady_min = nothing, qsteady_max = nothing)
    if beta === nothing 
        beta2 = QT(beta2_default(alg))
        beta1 = QT(beta1_default(alg, beta2))
        beta3 = QT(zero(beta1))
    else
        beta1, beta2, beta3 = beta
    end
    beta = SVector(map(float, promote(beta1, beta2, beta3))...)
    NewPIDController{QT, typeof(limiter)}(
        beta,
        QT(accept_safety),
        limiter,
        QT(qmin === nothing ? qmin_default(alg) : qmin),
        QT(qmax === nothing ? qmax_default(alg) : qmax),
        QT(qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min),
        QT(qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max),
    )
end

function Base.show(io::IO, controller::NewPIDController)
    print(io, "NewPIDController(beta=", controller.beta,
        ", accept_safety=", controller.accept_safety,
        ", limiter=", controller.limiter,
        ")")
end

mutable struct PIDControllerCache{T, Limiter} <: AbstractControllerCache
    controller::NewPIDController{T, Limiter}
    err::MVector{3, T} # history of the error estimates
    dt_factor::T
end

function setup_controller_cache(alg_cache, controller::NewPIDController{QT}) where QT
    err = MVector{3, QT}(true, true, true)
    PIDControllerCache(
        controller,
        err,
        QT(1 // 10^4),
    )
end

@inline function stepsize_controller!(integrator, cache::PIDControllerCache, alg)
    @unpack controller = cache
    beta1, beta2, beta3 = controller.beta
    @unpack qmax = controller

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
    @unpack controller = cache
    @unpack qsteady_min, qsteady_max = controller

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
    integrator.dt *= cache.qold
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
struct PredictiveController <: AbstractController
end

@inline function stepsize_controller!(integrator, controller::PredictiveController, alg)
    @unpack qmin, qmax, gamma = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)
    if iszero(EEst)
        q = inv(qmax)
    else
        if fac_default_gamma(alg)
            fac = gamma
        else
            if isfirk(alg)
                @unpack iter = integrator.cache
                @unpack maxiters = alg
            else
                @unpack iter, maxiters = integrator.cache.nlsolver
            end
            fac = min(gamma, (1 + 2 * maxiters) * gamma / (iter + 2 * maxiters))
        end
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = fastpower(EEst, expo) / fac
        @fastmath q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        integrator.qold = q
    end
    q
end

function step_accept_controller!(integrator, controller::PredictiveController, alg, q)
    @unpack qmin, qmax, gamma, qsteady_min, qsteady_max = integrator.opts

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
    integrator.erracc = max(1e-2, EEst)

    return integrator.dt / qacc
end

function step_reject_controller!(integrator, controller::PredictiveController, alg)
    @unpack dt, success_iter, qold = integrator
    integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
end


struct NewPredictiveController{T} <: AbstractController
    qmin::T
    qmax::T
    gamma::T
    qsteady_min::T
    qsteady_max::T
end

function NewPredictiveController(alg; kwargs...)
    PIController(Float64, alg; kwargs...)
end

function NewPredictiveController(QT, alg; qmin = nothing, qmax = nothing, gamma = nothing, qsteady_min = nothing, qsteady_max = nothing)
    NewPredictiveController{QT}(
        qmin === nothing ? qmin_default(alg) : qmin,
        qmax === nothing ? qmax_default(alg) : qmax,
        gamma === nothing ? gamma_default(alg) : gamma,
        qsteady_min === nothing ? qsteady_min_default(alg) : qsteady_min,
        qsteady_max === nothing ? qsteady_max_default(alg) : qsteady_max,
    )
end

mutable struct PredictiveControllerCache{T} <: AbstractControllerCache
    controller::NewPredictiveController{T}
    dtacc::T
    erracc::T
    qold::T
end

function setup_controller_cache(alg_cache, controller::NewPredictiveController{T}) where T
    PredictiveControllerCache{T}(
        controller,
        T(1),
        T(1),
        T(1),
    )
end

@inline function stepsize_controller!(integrator, cache::PredictiveControllerCache, alg)
    @unpack qmin, qmax, gamma = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)
    if iszero(EEst)
        q = inv(qmax)
    else
        if fac_default_gamma(alg)
            fac = gamma
        else
            if isfirk(alg)
                @unpack iter = integrator.cache
                @unpack maxiters = alg
            else
                @unpack iter, maxiters = integrator.cache.nlsolver
            end
            fac = min(gamma, (1 + 2 * maxiters) * gamma / (iter + 2 * maxiters))
        end
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = fastpower(EEst, expo) / fac
        @fastmath q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        cache.qold = q
    end
    cache.q = q
    q
end

function step_accept_controller!(integrator, cache::PredictiveControllerCache, alg, q)
    @assert q ≈ cache.q "Controller cache went out of sync with time stepping logic."
    @unpack dtacc, erracc, controller = cache
    @unpack qmin, qmax, gamma, qsteady_min, qsteady_max = controller

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
    cache.dtacc = integrator.dt
    cache.erracc = max(1e-2, EEst)

    return integrator.dt / qacc
end

function step_reject_controller!(integrator, cache::PredictiveControllerCache, alg)
    @unpack dt, success_iter = integrator
    @unpack qold = cache
    integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
end

# For a composite algorithm the default strategy is to switch forth and back between the controllers of the individual algorithms
struct CompositeController{T} <: AbstractController
    controllers::T
end

struct CompositeControllerCache{T} <: AbstractControllerCache
    caches::T
end

function setup_controller_cache(alg_cache::CompositeCache, cc::CompositeController)
    CompositeControllerCache(
        map((cache,controller)->default_controller_v7(cache, controller), alg_cache.caches, cc.controllers)
    )
end

@inline function accept_step_controller(integrator, cache::CompositeControllerCache)
    current_idx = integrator.cache.current
    accept_step_controller(integrator, cache.caches[current_idx])
end

@inline function stepsize_controller!(integrator, cache::CompositeControllerCache)
    current_idx = integrator.cache.current
    stepsize_controller!(integrator, cache.caches[current_idx])
end

@inline function step_accept_controller!(integrator, cache::CompositeControllerCache, q)
    current_idx = integrator.cache.current
    step_accept_controller!(integrator, cache.caches[current_idx], q)
end

@inline function step_reject_controller!(integrator, cache::CompositeControllerCache)
    current_idx = integrator.cache.current
    step_reject_controller!(integrator, cache.caches[current_idx])
end
