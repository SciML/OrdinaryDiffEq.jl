using StochasticDiffEq, DiffEqDevTools, Test
using Random
using SDEProblemLibrary: prob_sde_linear

f_linear_iip(du, u, p, t) = @.(du = 1.01 * u)
σ_linear_iip(du, u, p, t) = @.(du = 0.87 * u)

linear_analytic(u0, p, t, W) = @.(u0 * exp(0.63155t + 0.87W))

prob_sde_linear_iip = SDEProblem(
    SDEFunction(
        f_linear_iip, σ_linear_iip,
        analytic = linear_analytic
    ), [1 / 2], (0.0, 1.0)
)

Random.seed!(100)
dts = 1 .// 2 .^ (7:-1:3) #14->7 good plot

println("SROCKC2")
@time sim = test_convergence(
    dts, prob_sde_linear_iip, SROCKC2(), save_everystep = false, trajectories = Int(3.0e6),
    weak_timeseries_errors = false
)
@show sim.𝒪est[:weak_final]
@test abs(sim.𝒪est[:weak_final] - 2) < 0.4
#@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
sim = nothing

prob = prob_sde_linear

println("SROCKC2")
@time sim = test_convergence(
    dts, prob, SROCKC2(), save_everystep = false, trajectories = Int(1.0e6)
)
@show sim.𝒪est[:weak_final]
@test abs(sim.𝒪est[:weak_final] - 2) < 0.4
#@test abs(sim.𝒪est[:weak_l2]-2) < 0.3
#@test abs(sim.𝒪est[:weak_l∞]-2) < 0.3
