"""
 Tests for https://arxiv.org/abs/1303.5103 with test problems as in the paper.
 DRI1, RI1, RI3, RI5, RI6, RDI1WM, RDI2WM, RDI3WM, RDI4WM
"""

import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
using DiffEqDevTools
#using DiffEqGPU

function prob_func(prob, i, repeat)
    return remake(prob, seed = seeds[i])
end

"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"

u₀ = [1.0, 1.0]
function f2!(du, u, p, t)
    du[1] = -273 // 512 * u[1]
    return du[2] = -1 // 160 * u[1] - (-785 // 512 + sqrt(2) / 8) * u[2]
end
function g2!(du, u, p, t)
    du[1, 1] = 1 // 4 * u[1]
    du[1, 2] = 1 // 16 * u[1]
    du[2, 1] = (1 - 2 * sqrt(2)) / 4 * u[1]
    return du[2, 2] = 1 // 10 * u[1] + 1 // 16 * u[2]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 3.0)

h2(z) = z^2 # but apply it only to u[1]

prob = SDEProblem(f2!, g2!, u₀, tspan, noise_rate_prototype = zeros(2, 2))
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h2(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(1.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, DRI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("DRI1:", sim.𝒪est[:weak_final])

numtraj = Int(2.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
sim = test_convergence(
    dts, ensemble_prob, RI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order 3.43
println("RI1:", sim.𝒪est[:weak_final])

numtraj = Int(2.0e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
sim = test_convergence(
    dts, ensemble_prob, RI3(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order 2.46
println("RI3:", sim.𝒪est[:weak_final])

numtraj = Int(2.0e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
sim = test_convergence(
    dts, ensemble_prob, RI5(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order 2.57
println("RI5:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
sim = test_convergence(
    dts, ensemble_prob, RI6(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI6:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e5)
seed = 55
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI1WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test_broken abs(sim.𝒪est[:weak_final] - 2.0) < 0.3 # seems closer to 1.5?
@test abs(sim.𝒪est[:weak_final] - 1.5) < 0.3 # seems closer to 1.5?
println("RDI1WM:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI2WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 2.0) < 0.3
println("RDI2WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI3WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 3) < 0.4
println("RDI3WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI4WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 3) < 0.4
println("RDI4WM:", sim.𝒪est[:weak_final])
