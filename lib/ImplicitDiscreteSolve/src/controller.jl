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
Base.@kwdef struct KantorovichTypeController <: OrdinaryDiffEqCore.AbstractController
    Θmin::Float64
    p::Int64
    Θreject::Float64 = 0.95
    Θbar::Float64 = 0.5
    γ::Float64 = 0.95
    qmin::Float64 = 1 / 5
    qmax::Float64 = 5.0
    strict::Bool = true
end

function OrdinaryDiffEqCore.default_controller(
        alg::IDSolve, cache::IDSolveCache, _1, _2, _3
    )
    return KantorovichTypeController(; Θmin = 1 // 8, p = 1)
end

function OrdinaryDiffEqCore.stepsize_controller!(
        integrator, controller::KantorovichTypeController, alg::IDSolve
    )
    @inline g(x) = √(1 + 4x) - 1

    # Adapt dt with a priori estimate (Eq. 5.24)
    (; Θks) = integrator.cache
    (; Θbar, γ, Θmin, qmin, qmax, p) = controller

    Θ₀ = length(Θks) > 0 ? max(first(Θks), Θmin) : Θmin
    q = clamp(γ * (g(Θbar) / (g(Θ₀)))^(1 / p), qmin, qmax)

    return q
end

function OrdinaryDiffEqCore.step_accept_controller!(
        integrator, controller::KantorovichTypeController, alg::IDSolve, q
    )
    return q * integrator.dt
end

function OrdinaryDiffEqCore.step_reject_controller!(
        integrator, controller::KantorovichTypeController, alg::IDSolve
    )
    @inline g(x) = √(1 + 4x) - 1

    # Shorten dt according to (Eq. 5.24)
    (; Θks) = integrator.cache
    (; Θbar, Θreject, γ, Θmin, qmin, qmax, p) = controller
    for Θk in Θks
        if Θk > Θreject
            q = clamp(γ * (g(Θbar) / g(Θk))^(1 / p), qmin, qmax)
            integrator.dt = q * integrator.dt
            return
        end
    end
    return
end

function OrdinaryDiffEqCore.accept_step_controller(integrator, controller::KantorovichTypeController)
    (; Θks) = integrator.cache
    if controller.strict
        return all(controller.Θreject .< Θks)
    else
        return true
    end
end
