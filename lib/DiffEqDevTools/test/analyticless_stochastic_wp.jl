using Random
using StochasticDiffEq, DiffEqDevTools, Test
using SDEProblemLibrary: prob_sde_additivesystem

prob = prob_sde_additivesystem
prob = SDEProblem(prob.f, prob.g, prob.u0, (0.0, 0.1), prob.p)

reltols = 1.0 ./ 10.0 .^ (1:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [Dict(:alg => SRIW1())
          Dict(:alg => EM(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => RKMil(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => SRIW1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => SRA1(), :dts => 1.0 ./ 5.0 .^ ((1:length(reltols)) .+ 1),
              :adaptive => false)
          Dict(:alg => SRA1())]
_names = ["SRIW1", "EM", "RKMil", "SRIW1 Fixed", "SRA1 Fixed", "SRA1"]
test_dt = 0.1
wp = WorkPrecisionSet(prob, abstols, reltols, setups, test_dt;
    numruns = 5, names = _names, error_estimate = :l2)

se = get_sample_errors(prob, setups[1], numruns = 100, solution_runs = 100)
se = get_sample_errors(prob, setups[1], numruns = [5, 10, 25, 50, 100, 1000],
    solution_runs = 100)

println("Now weak error without analytical solution")

prob2 = SDEProblem((du, u, p, t) -> prob.f(du, u, p, t), prob.g, prob.u0, (0.0, 0.1),
    prob.p)
test_dt = 1 / 10^4
appxsol_setup = Dict(:alg => SRIW1(), :abstol => 1e-4, :reltol => 1e-4)
wp = WorkPrecisionSet(prob2, abstols, reltols, setups, test_dt;
    appxsol_setup = appxsol_setup,
    numruns = 5, names = _names, error_estimate = :weak_final)

println("Get sample errors")

se2 = get_sample_errors(prob2, setups[1], test_dt, appxsol_setup = appxsol_setup,
    numruns = [5, 10, 25, 50, 100], solution_runs = 20)

@test all(se[1:5] - se2 .< 1e-1)

# Ensemble Problem with non-commutative noise process

function prob_func(prob, i, repeat)
    remake(prob, seed = seeds[i])
end

u₀ = [1.0, 1.0]
function f2!(du, u, p, t)
    du[1] = -273 // 512 * u[1]
    du[2] = -1 // 160 * u[1] - (-785 // 512 + sqrt(2) / 8) * u[2]
    return nothing
end
function g2!(du, u, p, t)
    du[1, 1] = 1 // 4 * u[1]
    du[1, 2] = 1 // 16 * u[1]
    du[2, 1] = (1 - 2 * sqrt(2)) / 4 * u[1]
    du[2, 2] = 1 // 10 * u[1] + 1 // 16 * u[2]
    return nothing
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

reltols = 1.0 ./ 4.0 .^ (1:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [
    Dict(:alg => EM(), :dts => dts),
    Dict(:alg => SimplifiedEM(), :dts => dts),
    Dict(:alg => DRI1(), :dts => dts, :adaptive => false)
]
test_dt = 1 // 1000
appxsol_setup = Dict(:alg => EM(), :dt => test_dt)

# without analytical expectation value
wp1 = @time WorkPrecisionSet(ensemble_prob, abstols, reltols, setups, test_dt;
    maxiters = 1e7,
    verbose = false, save_everystep = false, save_start = false,
    appxsol_setup = appxsol_setup,
    trajectories = numtraj, error_estimate = :weak_final)

# with analytical expectation value
wp2 = @time WorkPrecisionSet(ensemble_prob, abstols, reltols, setups, test_dt;
    maxiters = 1e7,
    verbose = false, save_everystep = false, save_start = false,
    appxsol_setup = appxsol_setup, expected_value = exp(-3.0),
    trajectories = numtraj, error_estimate = :weak_final)

err1 = [wp1.wps[i].errors for i in 1:length(setups)]
err2 = [wp2.wps[i].errors for i in 1:length(setups)]

@test isapprox(err1, err2, atol = 1e-3)
