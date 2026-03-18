"""
 Tests for IRI1 (Implicit Rößler 1) drift-implicit weak order 2 SRK method.
 Based on test problems from https://arxiv.org/abs/1303.5103
"""

using StochasticDiffEq, DiffEqDevTools, Test
using Random

seed = 103478
function prob_func(prob, i, repeat)
    return remake(prob, seed = seeds[i])
end

@info "IRI1: Scalar oop weak convergence test"

# Test problem from the paper (scalar, oop)
numtraj = Int(2.0e6)
u₀ = 0.0
f(u, p, t) = 1 // 2 * u + sqrt(u^2 + 1)
g(u, p, t) = sqrt(u^2 + 1)
# Use dt values [1/16, 1/8, 1/4] - coarser values that show convergence
# Note: dt=1/2 causes nlsolver convergence failure for this nonlinear problem with theta=1
dts = 1 .// 2 .^ (4:-1:2)
tspan = (0.0, 2.0)

h1(z) = z^3 - 6 * z^2 + 8 * z
# analytical_sol(t) = E(h1(arsinh(X(t)))) = t^3-3*t^2+2*t
# analytical_sol(2) = 0

Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f, g, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(asinh(sol.u[end])), false),
    prob_func = prob_func
)

println("IRI1 weak convergence test")
sim = test_convergence(
    dts, ensemble_prob, IRI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@show sim.𝒪est[:weak_final]
# IRI1 has weak order 2
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("IRI1:", sim.𝒪est[:weak_final])
