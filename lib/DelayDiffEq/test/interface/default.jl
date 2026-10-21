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

stiff_dde(u, h, p, t) = -1000.0 * (u - h(p, t - 1.0))
function stiff_dde!(du, u, h, p, t)
    du[1] = -1000.0 * (u[1] - h(p, t - 1.0)[1])
    return nothing
end
stiff_probs = (
    DDEProblem(stiff_dde, 1.0, (p, t) -> 0.5, (0.0, 5.0); constant_lags = [1.0]),
    DDEProblem(stiff_dde!, [1.0], (p, t) -> [0.5], (0.0, 5.0); constant_lags = [1.0]),
)
@testset "Switch to $stiffalg, iip = $(isinplace(prob))" for prob in stiff_probs,
        (reltol, stiffalg) in (
            (1.0e-3, OrdinaryDiffEqDefault.DefaultSolverChoice.Rosenbrock23),
            (1.0e-8, OrdinaryDiffEqDefault.DefaultSolverChoice.Rodas5P),
        )

    sol = solve(prob; reltol)
    @test sol.retcode == ReturnCode.Success
    @test Int(stiffalg) in sol.alg_choice
    @test first(sol.u[end]) ≈ 0.5 atol = 1.0e-6
end
