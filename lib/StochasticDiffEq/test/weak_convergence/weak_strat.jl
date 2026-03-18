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

function prob_func(prob, i, repeat)
    return remake(prob, seed = seeds[i])
end

"""
 Test Scalar SDEs (oop)
"""

@info "Scalar oop noise"

numtraj = Int(8.0e6) # in the paper they use 1e9
u₀ = 0.1
p = [1.5, 0.1]
f(u, p, t) = p[1] * u
g(u, p, t) = p[2] * u
dts = 1 .// 2 .^ (5:-1:0)
tspan = (0.0, 1.0)

h1(z) = z

seed = 250
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f, g, u₀, tspan, p)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(sol.u[end]), false),
    prob_func = prob_func
)

sim = test_convergence(
    dts, ensemble_prob, RS1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("RS1:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, RS2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.39 #order is 2.36
println("RS2:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, NON(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.5
println("NON:", sim.𝒪est[:weak_final])

numtraj = Int(5.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
sim = test_convergence(
    dts, ensemble_prob, NON2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("NON2:", sim.𝒪est[:weak_final])

dts = 1 .// 2 .^ (4:-1:0)
numtraj = Int(1.0e6)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, COM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("COM:", sim.𝒪est[:weak_final])

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

f!(du, u, p, t) = du[1] = p[1] * u[1]
g!(du, u, p, t) = du[1] = p[2] * u[1]

numtraj = Int(5.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
dts = 1 .// 2 .^ (5:-1:0)

prob = SDEProblem(f!, g!, [u₀], tspan, p)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(sol.u[end][1]), false),
    prob_func = prob_func
)

sim = test_convergence(
    dts, ensemble_prob, RS1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.4
println("RS1:", sim.𝒪est[:weak_final])

numtraj = Int(2.0e7)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RS2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.45 #order is 2.415
println("RS2:", sim.𝒪est[:weak_final])

numtraj = Int(5.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, NON(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("NON:", sim.𝒪est[:weak_final])

numtraj = Int(5.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
sim = test_convergence(
    dts, ensemble_prob, NON2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("NON2:", sim.𝒪est[:weak_final])

dts = 1 .// 2 .^ (4:-1:0)
numtraj = Int(1.0e6)
seed = 10
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, COM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ * exp(1.0 * (p[1] + 0.5 * p[2]^2))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("COM:", sim.𝒪est[:weak_final])

"""
 Test Diagonal noise SDEs (iip)
"""

@info "Diagonal noise"

u₀ = [0.1, 0.1]
function f3!(du, u, p, t)
    du[1] = 299 // 200 * u[1]
    return du[2] = 299 // 200 * u[2]
end
function g3!(du, u, p, t)
    du[1] = 1 // 10 * u[1]
    return du[2] = 1 // 10 * u[2]
end
dts = 1 .// 2 .^ (4:-1:0)
tspan = (0.0, 1.0)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f3!, g3!, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h3(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(5.0e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RS1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.4 # order is 1.67
println("RS1:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e5)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, RS2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("RS2:", sim.𝒪est[:weak_final])

numtraj = Int(5.0e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, NON(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("NON:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, NON2(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35
println("NON2:", sim.𝒪est[:weak_final])

numtraj = Int(1.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, COM(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("COM:", sim.𝒪est[:weak_final])
