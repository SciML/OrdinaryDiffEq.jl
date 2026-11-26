# JVODE

"""
    JVODEController(;
        qmin, qmax, qsteady_min, qsteady_max, gamma,
        qmax_first_step, failfactor
    )

Step-size controller for the variable-order Nordsieck-form `JVODE` family.
Composes the standard step-size knobs via [`CommonControllerOptions`](@ref); the
adaptive logic is integrated into the algorithm itself, so the controller
cache falls through to alg-level dispatch but the knobs are exposed as a
real, settable controller (instead of being hard-wired on the algorithm
struct as before).
"""
struct JVODEController{B <: Union{NamedTuple, CommonControllerOptions}} <: AbstractController
    basic::B
end

JVODEController(; kwargs...) = JVODEController(NamedTuple(kwargs))

JVODEController(alg; kwargs...) = JVODEController(Float64, alg; kwargs...)
JVODEController(::Type{QT}, alg; kwargs...) where {QT} =
    JVODEController(resolve_basic(NamedTuple(kwargs), alg, QT))

mutable struct JVODEControllerCache{J <: JVODEController, C, E} <: AbstractControllerCache
    controller::J
    cache::C
    EEst::E
end

function setup_controller_cache(
        alg::JVODE, cache, controller::JVODEController, ::Type{E}, disco_probs,
    ) where {E}
    QT = _resolved_QT(controller.basic)
    resolved = JVODEController(resolve_basic(controller.basic, alg, QT; disco_probs))
    return JVODEControllerCache(resolved, cache, oneunit(E))
end

function setup_controller_cache(
        alg::JVODE, cache, ::DummyController, ::Type{E}, disco_probs,
    ) where {E}
    basic = (;
        qmin = alg.qmin, qmax = alg.qmax,
        qsteady_min = alg.qsteady_min, qsteady_max = alg.qsteady_max,
        qmax_first_step = qmax_first_step_default(alg), gamma = gamma_default(alg),
        failfactor = failfactor_default(alg), discontinuity_detection = false, disco_probs,
    )
    return JVODEControllerCache(JVODEController(basic), cache, oneunit(E))
end

# Algorithm owns the stepsize logic; controller cache delegates back to
# alg-level dispatch (mirroring how DummyControllerCache used to behave).
@inline stepsize_controller!(integrator, ::JVODEControllerCache, alg) =
    stepsize_controller!(integrator, alg)
@inline step_accept_controller!(integrator, ::JVODEControllerCache, alg, q) =
    step_accept_controller!(integrator, alg, q)
@inline step_reject_controller!(integrator, ::JVODEControllerCache, alg) =
    step_reject_controller!(integrator, alg)
@inline accept_step_controller(
    integrator, cache::JVODEControllerCache, alg,
) = get_EEst(cache) <= 1

function stepsize_controller!(integrator, alg::JVODE)
    if iszero(get_EEst(integrator))
        η = get_qmax(integrator)
    else
        η = integrator.cache.η
        integrator.cache.ηold = η
    end
    return η
end

function step_accept_controller!(integrator, alg::JVODE, η)
    q = inv(η)
    if q <= get_qsteady_max(integrator) && q >= get_qsteady_min(integrator)
        q = one(q)
    end
    return integrator.dt / q  # dtnew
end

function step_reject_controller!(integrator, alg::JVODE)
    return integrator.dt *= integrator.cache.η
end
