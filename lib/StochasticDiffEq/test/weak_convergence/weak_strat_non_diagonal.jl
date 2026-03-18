"""
 Tests for  https://link.springer.com/article/10.1007/s10543-007-0130-3 with test problems as in the paper.
 RS1, RS2
 and for https://www.sciencedirect.com/science/article/pii/S0377042706003906
 NON, COM, NON2
"""

import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
using DiffEqDevTools
#using DiffEqGPU

seed = 137475

function prob_func(prob, i, repeat)
    return remake(prob, seed = seeds[i])
end

"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"

u₀ = [0.1, 0.1]
function f2!(du, u, p, t)
    du[1] = 5 // 4 * u[2] - 5 // 4 * u[1]
    return du[2] = 1 // 4 * u[1] - 1 // 4 * u[2]
end
function g2!(du, u, p, t)
    du[1, 1] = sqrt(3) / 2 * (u[1] - u[2])
    du[1, 2] = 1 // 2 * (u[1] + u[2])
    #du[2,1] = 0
    return du[2, 2] = u[1]
end
dts = 1 .// 2 .^ (4:-1:1)
tspan = (0.0, 1.0)

h2(z) = z^2 # E(x_i) = 1/10 exp(1/2t) or E(x_1* x_2) = 1/100 exp(2t)

prob = SDEProblem(f2!, g2!, u₀, tspan, noise_rate_prototype = zeros(2, 2))
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h2(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(1.0e8)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RS1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(2)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RS1:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RS2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(2)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RS2:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, NON(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(2)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("NON:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, NON2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(2)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("NON2:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e8)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, COM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(2)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("COM:", sim.𝒪est[:weak_final])
# COM tests are passing; problem might be not hard enough..
