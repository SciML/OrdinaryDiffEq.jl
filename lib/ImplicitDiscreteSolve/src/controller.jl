"""
    KantorovichTypeController

Default controller for implicit discrete solvers. Assuming a Newton method is used
to solve the nonlinear problem, this controller uses convergence rate estimates to
adapt the step size based on an posteriori estimate of the change in the Newton
convergence radius between steps as described below.

Given the convergence rate estimate Θ₀ for the first iteration, the step size controller 
adapts the time step as `dtₙ₊₁ = γ * (g(Θbar) / (g(Θ₀)))^(1 / p) dtₙ`
with `g(x) = √(1 + 4x) - 1`. `p` denotes the order of the solver -- i.e. the order of
the extrapolation algorithm to compute the initial guess for the solve at `tₙ` given a 
solution at `tₙ₋₁` -- and `Θbar` denotes the  target convergence rate. `γ` is a safety factor.

A factor `Θreject` controls the rejection of a time step if any `Θₖ > Θreject`.
In this case the first Θₖ violating this criterion is taken and the time step is adapted
such that `dtₙ₊₁ = γ * (g(Θbar) / (g(Θk)))^(1 / p) dtₙ`. This behavior can be changed
by setting `strict = false`. In this case the step is accepted whenever the Newton
solver converges.

The controller furthermore limits the growth and shrinkage of the time step by a factor
between `qmin` and `qmax`.

The baseline algorithm has been derived in Peter Deuflhard's book "Newton Methods for
Nonlinear Problems" in Section 5.1.3 (Adaptive pathfollowing algorithms). Please note
that some implementation details deviate from the original algorithm.
"""
struct KantorovichTypeController{B <: Union{NamedTuple, OrdinaryDiffEqCore.CommonControllerOptions}, T} <: AbstractController
    basic::B
    Θmin::T
    p::Int64
    Θreject::T
    Θbar::T
    γ::T
    strict::Bool
end

function KantorovichTypeController(;
        Θmin, p, Θreject = 0.95, Θbar = 0.5, γ = 0.95,
        qmin = 1 // 5, qmax = 5, strict = true,
        kwargs...,
    )
    T = promote_type(typeof(Θmin), typeof(Θreject), typeof(Θbar), typeof(γ))
    basic = (; qmin, qmax, kwargs...)
    return KantorovichTypeController{typeof(basic), T}(
        basic, T(Θmin), Int64(p), T(Θreject), T(Θbar), T(γ), strict,
    )
end

mutable struct KantorovichTypeControllerCache{T, E, NLPType} <: AbstractControllerCache
    controller::KantorovichTypeController{OrdinaryDiffEqCore.CommonControllerOptions{T, NLPType}, T}
    EEst::E
end

function OrdinaryDiffEqCore.default_controller(
        QT, alg::IDSolve,
    )
    return KantorovichTypeController(; Θmin = QT(1 // 8), p = 1)
end

function OrdinaryDiffEqCore.setup_controller_cache(alg, cache, controller::KantorovichTypeController, ::Type{E}, disco_probs) where {E}
    QT = OrdinaryDiffEqCore._resolved_QT(controller.basic)
    basic = OrdinaryDiffEqCore.resolve_basic(controller.basic, alg, QT; disco_probs)
    resolved = KantorovichTypeController{typeof(basic), QT}(
        basic, QT(controller.Θmin), controller.p,
        QT(controller.Θreject), QT(controller.Θbar), QT(controller.γ), controller.strict,
    )
    T = QT
    return KantorovichTypeControllerCache{T, E, eltype(disco_probs)}(resolved, oneunit(E))
end

function OrdinaryDiffEqCore.stepsize_controller!(
        integrator, cache::KantorovichTypeControllerCache, alg::IDSolve
    )
    (; controller) = cache
    @inline g(x) = √(1 + 4x) - 1

    # Adapt dt with a priori estimate (Eq. 5.24)
    (; Θks) = integrator.cache
    (; Θbar, γ, Θmin, p) = controller
    (; qmin, qmax) = controller.basic

    Θ₀ = length(Θks) > 0 ? max(first(Θks), Θmin) : Θmin
    q = clamp(γ * (g(Θbar) / (g(Θ₀)))^(1 / p), qmin, qmax)

    return q
end

function OrdinaryDiffEqCore.step_accept_controller!(
        integrator, cache::KantorovichTypeControllerCache, alg::IDSolve, q
    )
    return q * integrator.dt
end

function OrdinaryDiffEqCore.step_reject_controller!(
        integrator, cache::KantorovichTypeControllerCache, alg::IDSolve
    )
    (; controller) = cache
    @inline g(x) = √(1 + 4x) - 1

    # Shorten dt according to (Eq. 5.24)
    (; Θks) = integrator.cache
    (; Θbar, Θreject, γ, Θmin, p) = controller
    (; qmin, qmax) = controller.basic
    for Θk in Θks
        if Θk > Θreject
            q = clamp(γ * (g(Θbar) / g(Θk))^(1 / p), qmin, qmax)
            integrator.dt = q * integrator.dt
            return
        end
    end
    return
end

function _accept_kantorovich_step(controller, Θks)
    return !controller.strict || all(Θk -> Θk <= controller.Θreject, Θks)
end

function OrdinaryDiffEqCore.accept_step_controller(integrator, cache::KantorovichTypeControllerCache, alg)
    (; controller) = cache
    (; Θks) = integrator.cache
    return _accept_kantorovich_step(controller, Θks)
end

function OrdinaryDiffEqCore.sync_controllers!(cache1::KantorovichTypeControllerCache, cache2::KantorovichTypeControllerCache)
    return nothing
end
