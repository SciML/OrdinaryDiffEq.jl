"""
 Tests for SIE and SME methods
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

# PC exercise 14.2.2
numtraj = Int(1.0e5)
u₀ = 0.1
f(u, p, t) = p[1] * u
g(u, p, t) = p[2] * u
dts = 1 .// 2 .^ (6:-1:3)
tspan = (0.0, 1.0)
p = [3 // 2, 1 // 100]

h1(z) = z
#analytical_sol(t) = E(f(X(t))) =

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f, g, u₀, tspan, p)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(sol.u[end]), false),
    prob_func = prob_func
)

sim = test_convergence(
    dts, ensemble_prob, SIEA(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.35
println("SIEA:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SMEA(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.34
println("SMEA:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SIEB(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.35
println("SIEB:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SMEB(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.34
println("SMEB:", sim.𝒪est[:weak_final])

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

u₀ = [0.1]
f1!(du, u, p, t) = (du[1] = p[1] * u[1])
g1!(du, u, p, t) = (du[1] = p[2] * u[1])
dts = 1 .// 2 .^ (6:-1:3)
tspan = (0.0, 1.0)
p = [3 // 2, 1 // 100]

prob = SDEProblem(f1!, g1!, u₀, tspan, p)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h1(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(1.0e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, SIEA(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.35
println("SIEA:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SMEA(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.34
println("SMEA:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SIEB(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.35
println("SIEB:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SMEB(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = u₀ .* exp(1.0 * (p[1]))
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.35 # order is 2.34
println("SMEB:", sim.𝒪est[:weak_final])

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
dts = 1 .// 2 .^ (5:-1:1)
tspan = (0.0, 1.0)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f3!, g3!, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h3(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(1.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, SIEA(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("SIEA:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SMEA(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("SMEA:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SIEB(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("SIEB:", sim.𝒪est[:weak_final])

sim = test_convergence(
    dts, ensemble_prob, SMEB(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("SMEB:", sim.𝒪est[:weak_final])
