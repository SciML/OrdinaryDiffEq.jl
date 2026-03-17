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
 Test Scalar SDEs (oop)
"""

@info "Scalar oop noise"

numtraj = Int(1.0e7) # in the paper they use 1e9
u₀ = 0.0
f(u, p, t) = 1 // 2 * u + sqrt(u^2 + 1)
g(u, p, t) = sqrt(u^2 + 1)
dts = 1 .// 2 .^ (4:-1:1)
tspan = (0.0, 2.0) # 2.0 in paper

h1(z) = z^3 - 6 * z^2 + 8 * z
#analytical_sol(t) = E(f(X(t))) = E(h1(arsinh(X(t))) = t^3-3*t^2+2*t
#analytical_sol(2) = 0 and analytical_sol(1)=0

Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f, g, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(asinh(sol.u[end])), false),
    prob_func = prob_func
)

sim = test_convergence(
    dts, ensemble_prob, DRI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("DRI1:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI1:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI3(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI3:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI5(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI5:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI6(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI6:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI1WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
println("RDI1WM:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI2WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RDI2WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI3WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RDI3WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI4WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.5 # order is 2.568 with 3e5 trajectories
println("RDI4WM:", sim.𝒪est[:weak_final])

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

u₀ = [0.0]
f1!(du, u, p, t) = @.(du = 1 // 2 * u + sqrt(u^2 + 1))
g1!(du, u, p, t) = @.(du = sqrt(u^2 + 1))
dts = 1 .// 2 .^ (4:-1:1)
tspan = (0.0, 2.0)

h1(z) = z^3 - 6 * z^2 + 8 * z

prob = SDEProblem(f1!, g1!, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(asinh(sol.u[end][1])), false),
    prob_func = prob_func
)

numtraj = Int(6.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, DRI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("DRI1:", sim.𝒪est[:weak_final])

numtraj = Int(8.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, DRI1NM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("DRI1NM:", sim.𝒪est[:weak_final])

numtraj = Int(9.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.4
println("RI1:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI3(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI3:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI5(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI5:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RI6(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RI6:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI1WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 1) < 0.3
println("RDI1WM:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e7)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RDI2WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("RDI2WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI3WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RDI3WM:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RDI4WM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3
println("RDI4WM:", sim.𝒪est[:weak_final])
