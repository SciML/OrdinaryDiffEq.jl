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
seed = 103473
function prob_func(prob, i, repeat)
    return remake(prob, seed = seeds[i])
end

"""
 Test Diagonal noise SDEs (iip), SIAM Journal on Numerical Analysis, 47 (2009), pp. 1713–1738
"""

@info "Diagonal noise"

u₀ = [0.1, 0.1]
function f3!(du, u, p, t)
    du[1] = 3 // 2 * u[1]
    return du[2] = 3 // 2 * u[2]
end
function g3!(du, u, p, t)
    du[1] = 1 // 10 * u[1]
    return du[2] = 1 // 10 * u[2]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 1.0)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f3!, g3!, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h3(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(5.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, DRI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 2.91
println("DRI1:", sim.𝒪est[:weak_final])

numtraj = Int(8.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, DRI1NM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.4 # order is 2.91
println("DRI1NM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 3.05
println("RI1:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RI3(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 2.77
println("RI3:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RI5(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.45
println("RI5:", sim.𝒪est[:weak_final])

numtraj = Int(3.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI6(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 1.5) < 0.5 # order is 1.55
println("R6:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI1WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.5 # order is 1.44
println("RDI1WM:", sim.𝒪est[:weak_final])

numtraj = Int(3.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI2WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 1.5) < 0.5 # order is 1.55
println("RDI2WM:", sim.𝒪est[:weak_final])

numtraj = Int(5.0e4)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI3WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 2.84
println("RDI3WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI4WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 2.91
println("RDI4WM:", sim.𝒪est[:weak_final])
