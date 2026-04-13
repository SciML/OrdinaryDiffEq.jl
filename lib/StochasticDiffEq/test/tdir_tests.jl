using StochasticDiffEq, Test, Random, DiffEqNoiseProcess
using SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear

Random.seed!(100)

for prob in (prob_sde_linear, prob_sde_2Dlinear)
    prob2 = remake(prob; tspan = (1.0, 0.0))

    solEM = solve(prob2, EM(), dt = -1 // 2^(4), tstops = [0.33])
    @test solEM.t[15] > solEM.t[16]

    solSRIW1 = solve(prob2, SRIW1(), tstops = [0.33])
    @test solSRIW1.t[15] > solSRIW1.t[16]
end

# Test reversing the noise

α = 1
β = 1
u₀ = 1 / 2
f(u, p, t) = α * u
g(u, p, t) = β * u
dt = 1 // 2^(4)
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u₀, (0.0, 1.0))
sol = solve(prob, EulerHeun(), dt = 0.01, save_noise = true)
_sol = deepcopy(sol) # to make sure the plot is correct
W3 = NoiseGrid(reverse!(_sol.t), reverse!(_sol.W))
prob3 = SDEProblem(f, g, sol.u[end], (1.0, 0.0), noise = W3)
sol2 = solve(prob3, EulerHeun(), dt = 0.01)
@test sol.u ≈ reverse!(sol2.u) atol = 1.0e-1
