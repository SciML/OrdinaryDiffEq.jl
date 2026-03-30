using OrdinaryDiffEq, StochasticDelayDiffEq, ParameterizedFunctions, Test, Random
using ParameterizedFunctions.ModelingToolkit # macro hygiene
f = @ode_def LotkaVolterra begin
    dx = 1.5x - x * y
    dy = -3y + x * y
end

prob = ODEProblem(f, big.([1.0; 1.0]), (big(0.0), big(10.0)))

using DiffEqDevTools

dts = big(1 / 2) .^ (8:-1:4)
test_setup = Dict(:alg => Vern9(), :reltol => 1e-25, :abstol => 1e-25)
sim1 = analyticless_test_convergence(dts, prob, Tsit5(), test_setup)
sim2 = analyticless_test_convergence(dts, prob, Vern9(), test_setup)

if VERSION >= v"1.8"
    @test_throws "No analytical solution set." test_convergence(dts, prob, Tsit5())
end

@test abs(sim1.𝒪est[:final] - 5) < 0.2
@test abs(sim2.𝒪est[:final] - 9) < 0.2

function f2(du, u, p, t)
    du .= 1.01u
end
function g2(du, u, p, t)
    du .= 1.01u
end
prob = SDEProblem(f2, g2, [1.0; 1.0], (0.0, 1.0))

using StochasticDiffEq

dts = (1 / 2) .^ (7:-1:3)
test_dt = 1 / 2^8
Random.seed!(100)
sim1 = analyticless_test_convergence(dts, prob, SRIW1(), test_dt, trajectories = 100)
@test abs(sim1.𝒪est[:final] - 1.5) < 0.4
@show sim1.𝒪est[:final]

dts = (1 / 2) .^ (7:-1:4)
test_dt = 1 / 2^8
sim2 = analyticless_test_convergence(dts, prob, SRIW1(), test_dt, trajectories = 100,
    use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.5) < 0.3
@show sim2.𝒪est[:final]

# EnsembleProblem

function prob_func(prob, i, repeat)
    remake(prob, seed = seeds[i])
end

u₀ = [1.0, 1.0]
function f2!(du, u, p, t)
    du[1] = -273 // 512 * u[1]
    du[2] = -1 // 160 * u[1] - (-785 // 512 + sqrt(2) / 8) * u[2]
end
function g2!(du, u, p, t)
    du[1, 1] = 1 // 4 * u[1]
    du[1, 2] = 1 // 16 * u[1]
    du[2, 1] = (1 - 2 * sqrt(2)) / 4 * u[1]
    du[2, 2] = 1 // 10 * u[1] + 1 // 16 * u[2]
end
dts = 1 .// 2 .^ (3:-1:0)
tspan = (0.0, 3.0)

h2(z) = z^2 # but apply it only to u[1]

prob = SDEProblem(f2!, g2!, u₀, tspan, noise_rate_prototype = zeros(2, 2))

numtraj = Int(1e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)
ensemble_prob = EnsembleProblem(prob;
    output_func = (sol, i) -> (h2(sol[1, end]), false),
    prob_func = prob_func)
sim = test_convergence(dts, ensemble_prob, DRI1(), save_everystep = false,
    trajectories = numtraj, save_start = false, adaptive = false,
    weak_timeseries_errors = false, weak_dense_errors = false,
    expected_value = exp(-3.0))

@test abs(sim.𝒪est[:weak_final] - 2.0) < 0.3
@show sim.𝒪est[:weak_final]

### SDDE

function hayes_modelf(du, u, h, p, t)
    τ, a, b, c, α, β, γ = p
    du .= a .* u .+ b .* h(p, t - τ) .+ c
end
function hayes_modelg(du, u, h, p, t)
    τ, a, b, c, α, β, γ = p
    du .= α .* u .+ β .* h(p, t - τ) .+ γ
end
h(p, t) = (ones(1) .+ t);
tspan = (0.0, 10.0)

pmul = [1.0, -4.0, -2.0, 10.0, -1.3, -1.2, 1.1]
padd = [1.0, -4.0, -2.0, 10.0, -0.0, -0.0, 0.1]

prob = SDDEProblem(hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
    constant_lags = (pmul[1],));
dts = (1 / 2) .^ (7:-1:3)
test_dt = 1 / 2^8
sim2 = analyticless_test_convergence(dts, prob, RKMil(), test_dt, trajectories = 100,
    use_noise_grid = false)
@test abs(sim2.𝒪est[:final] - 1.0) < 0.3
