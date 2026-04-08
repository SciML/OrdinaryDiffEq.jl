using DelayDiffEq
using StochasticDiffEqImplicit
using Test

# Hayes equation with implicit SDE methods (regression test for #3363)
begin
    function hayes_modelf(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        return du .= a .* u .+ b .* h(p, t - τ) .+ c
    end
    function hayes_modelg(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        return du .= α .* u .+ β .* h(p, t - τ) .+ γ
    end
    h(p, t) = (ones(1) .+ t)
    tspan = (0.0, 0.1)
end

pmul = [1.0, -4.0, -2.0, 10.0, -1.3, -1.2, 1.1]

prob = SDDEProblem(
    hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
    constant_lags = (pmul[1],)
)

# ImplicitEM with explicit dt
sol = solve(prob, MethodOfSteps(ImplicitEM()), dt = 0.01)
@test SciMLBase.successful_retcode(sol)
@test sol.u[end] != zeros(1)

# ImplicitEM without explicit dt (tests auto-dt with SDE algorithm order)
sol = solve(prob, MethodOfSteps(ImplicitEM()))
@test SciMLBase.successful_retcode(sol)
@test sol.u[end] != zeros(1)

# ImplicitEulerHeun
sol = solve(prob, MethodOfSteps(ImplicitEulerHeun()), dt = 0.01)
@test SciMLBase.successful_retcode(sol)
@test sol.u[end] != zeros(1)

# ImplicitRKMil
sol = solve(prob, MethodOfSteps(ImplicitRKMil()), dt = 0.01, adaptive = false)
@test SciMLBase.successful_retcode(sol)
@test sol.u[end] != zeros(1)

# ISSEM (Implicit Split-Step EM)
sol = solve(prob, MethodOfSteps(ISSEM()), dt = 0.01)
@test SciMLBase.successful_retcode(sol)
@test sol.u[end] != zeros(1)

# ISSEulerHeun (Implicit Split-Step Euler-Heun)
sol = solve(prob, MethodOfSteps(ISSEulerHeun()), dt = 0.01)
@test SciMLBase.successful_retcode(sol)
@test sol.u[end] != zeros(1)
