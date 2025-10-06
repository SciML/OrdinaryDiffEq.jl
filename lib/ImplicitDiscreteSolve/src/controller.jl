Base.@kwdef struct KantorovichTypeController <: OrdinaryDiffEqCore.AbstractController
    Θmin::Float64
    p::Int64
    Θreject::Float64 = 0.95
    Θbar::Float64 = 0.5
    γ::Float64    = 0.95
    qmin::Float64 = 1/5
    qmax::Float64 = 5.0
end

function OrdinaryDiffEqCore.default_controller(alg::IDSolve, cache::IDSolveCache, _1, _2, _3)
    return KantorovichTypeController(;Θmin=1//8, p=1)
end

function OrdinaryDiffEqCore.stepsize_controller!(
    integrator, controller::KantorovichTypeController, alg::IDSolve
)
    @inline g(x) = √(1+4x) - 1

    # Adapt dt with a priori estimate (Eq. 5.24)
    (; Θks) = integrator.cache
    (; Θbar, γ, Θmin, qmin, qmax, p) = controller

    Θ₀ = length(Θks) > 0 ? max(first(Θks), Θmin) : Θmin
    q = clamp(γ * (g(Θbar)/(g(Θ₀)))^(1/p), qmin, qmax)

    return q
end

function OrdinaryDiffEqCore.step_accept_controller!(
    integrator, controller::KantorovichTypeController, alg::IDSolve, q
)
    @info integrator.dt, q
    return q * integrator.dt
end

function OrdinaryDiffEqCore.step_reject_controller!(
    integrator, controller::KantorovichTypeController, alg::IDSolve
)
    @inline g(x) = √(1+4x) - 1

    # Shorten dt according to (Eq. 5.24)
    (; Θks) = cache.inner_solver_cache
    (; Θbar, Θreject, γ, Θmin, qmin, qmax, p) = controller
    for Θk in Θks
        if Θk > Θreject
            q = clamp(γ * (g(Θbar)/g(Θk))^(1/p), qmin, qmax)
            integrator.dt = q * integrator.dt
            return
        end
    end
end
