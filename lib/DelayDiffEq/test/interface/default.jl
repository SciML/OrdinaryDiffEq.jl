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

stiff_exact(t) = 0.5 + 0.25 * sin(t)
stiff_dde(u, h, p, t) = 0.25 * cos(t) +
    100.0 * (h(p, t - 1.0) - stiff_exact(t - 1.0)) - 1000.0 * (u - stiff_exact(t))
function stiff_dde!(du, u, h, p, t)
    du[1] = 0.25 * cos(t) +
        100.0 * (h(p, t - 1.0)[1] - stiff_exact(t - 1.0)) - 1000.0 * (u[1] - stiff_exact(t))
    return nothing
end
stiff_probs = (
    DDEProblem(stiff_dde, stiff_exact(0.0), (p, t) -> stiff_exact(t), (0.0, 5.0); constant_lags = [1.0]),
    DDEProblem(stiff_dde!, [stiff_exact(0.0)], (p, t) -> [stiff_exact(t)], (0.0, 5.0); constant_lags = [1.0]),
)
@testset "Switch to $stiffalg, iip = $(isinplace(prob))" for prob in stiff_probs,
        (reltol, stiffalg) in (
            (1.0e-3, OrdinaryDiffEqDefault.DefaultSolverChoice.Rosenbrock23),
            (1.0e-8, OrdinaryDiffEqDefault.DefaultSolverChoice.Rodas5P),
        )

    abstol = 1.0e-6
    sol = solve(prob; reltol, abstol)
    @test sol.retcode == ReturnCode.Success
    @test sol.t[end] == 5.0
    @test Int(stiffalg) in sol.alg_choice
    error(t) = abs(first(sol(t)) - stiff_exact(t))
    grid_error = maximum(error, 0.0:0.0005:5.0)
    interval_error = maximum(
        error(a + c * (b - a)) for (a, b) in zip(sol.t[1:(end - 1)], sol.t[2:end])
            for c in (0.25, 0.5, 0.75)
    )
    @test max(grid_error, interval_error) < 2 * (abstol + 0.75 * reltol)
end
