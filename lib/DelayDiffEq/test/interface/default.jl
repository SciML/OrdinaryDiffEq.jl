using DelayDiffEq, Test
using OrdinaryDiffEqDefault
function delay_lotka_volterra(du, u, h, p, t)
    # Model parameters.
    α, β, γ, δ = p
    # Current state.
    x, y = u
    # Evaluate differential equations
    du[1] = α * h(p, t - 1; idxs = 1) - β * x * y
    du[2] = -γ * y + δ * x * y
    return nothing
end
p = (1.5, 1.0, 3.0, 1.0)
u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
h(p, t; idxs::Int) = 1.0
prob_dde = DDEProblem(delay_lotka_volterra, u0, h, tspan, p);
alg = MethodOfSteps(DefaultODEAlgorithm())
sol_dde = solve(prob_dde, alg)
sol_dde2 = solve(prob_dde)

@test sol_dde == sol_dde2
