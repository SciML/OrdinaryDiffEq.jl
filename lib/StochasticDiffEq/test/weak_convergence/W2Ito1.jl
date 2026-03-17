"""
  Tests for W2Ito1 from Tang, X., & Xiao, A. (2017). Efficient weak second-order stochastic
  Runge–Kutta methods for Itô stochastic differential equations. BIT Numerical Mathematics,
  57, 241-260.
"""

import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
using DiffEqDevTools
seed = 103473
function prob_func(prob, i, repeat)
    return remake(prob, seed = seeds[i])
end

"""
 Test OOP
"""

@info "Scalar noise"

numtraj = Int(2.0e6) # in the paper they use 1e9
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
    dts, ensemble_prob, W2Ito1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("W2Ito1:", sim.𝒪est[:weak_final])

@info "Diagonal noise"

u₀ = [0.1, 0.1]
function f2(u, p, t)
    return [3 // 2 * u[1], 3 // 2 * u[2]]
end
function g2(u, p, t)
    return [1 // 10 * u[1], 1 // 10 * u[2]]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 1.0)

h2(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f2, g2, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h2(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(1.0e5)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, W2Ito1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 2.4
println("W2Ito1:", sim.𝒪est[:weak_final])

@info "Non-commutative noise"

u₀ = [1.0, 1.0]
function f3(u, p, t)
    return [-273 // 512 * u[1], -1 // 160 * u[1] - (-785 // 512 + sqrt(2) / 8) * u[2]]
end
function g3(u, p, t)
    return [
        1 // 4 * u[1] 1 // 16 * u[1]
        (1 - 2 * sqrt(2)) / 4 * u[1] 1 // 10 * u[1] + 1 // 16 * u[2]
    ]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 3.0)

h3(z) = z^2 # but apply it only to u[1]

prob = SDEProblem(f3, g3, u₀, tspan, noise_rate_prototype = zeros(2, 2))
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
    dts, ensemble_prob, W2Ito1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3

println("W2Ito1:", sim.𝒪est[:weak_final])

# """
#  Test IIP
# """

@info "Scalar noise"

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

numtraj = Int(1.0e6)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, W2Ito1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 0.0
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3
println("W2Ito1:", sim.𝒪est[:weak_final])

@info "Diagonal noise"

u₀ = [0.1, 0.1]
function f2!(du, u, p, t)
    du[1] = 3 // 2 * u[1]
    return du[2] = 3 // 2 * u[2]
end
function g2!(du, u, p, t)
    du[1] = 1 // 10 * u[1]
    return du[2] = 1 // 10 * u[2]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 1.0)

h2(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f2!, g2!, u₀, tspan)
ensemble_prob = EnsembleProblem(
    prob;
    output_func = (sol, i) -> (h2(sol.u[end][1]), false),
    prob_func = prob_func
)

numtraj = Int(1.0e5)
Random.seed!(seed)
seeds = rand(UInt, numtraj)

sim = test_convergence(
    dts, ensemble_prob, W2Ito1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = 1 // 100 * exp(301 // 100)
)
@test -(sim.𝒪est[:weak_final] - 2) < 0.3 # order is 2.4
println("W2Ito1:", sim.𝒪est[:weak_final])

@info "Non-commutative noise"

u₀ = [1.0, 1.0]
function f3!(du, u, p, t)
    du[1] = -273 // 512 * u[1]
    return du[2] = -1 // 160 * u[1] - (-785 // 512 + sqrt(2) / 8) * u[2]
end
function g3!(du, u, p, t)
    du[1, 1] = 1 // 4 * u[1]
    du[1, 2] = 1 // 16 * u[1]
    du[2, 1] = (1 - 2 * sqrt(2)) / 4 * u[1]
    return du[2, 2] = 1 // 10 * u[1] + 1 // 16 * u[2]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 3.0)

h3(z) = z^2 # but apply it only to u[1]

prob = SDEProblem(f3!, g3!, u₀, tspan, noise_rate_prototype = zeros(2, 2))
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
    dts, ensemble_prob, W2Ito1(),
    save_everystep = false, trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0)
)
@test abs(sim.𝒪est[:weak_final] - 2) < 0.3

println("W2Ito1:", sim.𝒪est[:weak_final])
