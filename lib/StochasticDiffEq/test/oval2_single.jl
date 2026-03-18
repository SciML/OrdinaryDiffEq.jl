using StochasticDiffEq, Random
using SDEProblemLibrary: oval2ModelExample
Random.seed!(200)
prob = oval2ModelExample(largeFluctuations = true, useBigs = false)
quick_prob = deepcopy(prob)
quick_prob.tspan = (0.0, 1.0)

using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.0001
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
@benchmark begin
    Random.seed!(100)
    sol = solve(
        quick_prob, SRIW1(), dt = (1 / 2)^(18), progress_steps = Int(1.0e5),
        adaptivealg = :RSwM3, progress = false, qmax = 4, save_everystep = false,
        saveat = 0.1, abstol = 1.0e-5, reltol = 1.0e-3
    )
end

using ProfileView
Random.seed!(100)
@profile sol = solve(
    quick_prob, SRIW1(), dt = (1 / 2)^(18), progress_steps = Int(1.0e5),
    adaptivealg = :RSwM3, progress = false, qmax = 4, save_everystep = false,
    saveat = 0.1, abstol = 1.0e-5, reltol = 1.0e-3
)

ProfileView.view()

Random.seed!(100)
@time sol = solve(
    quick_prob, SRIW1(), dt = (1 / 2)^(18), progress_steps = Int(1.0e5),
    progress = false, qmax = 4, save_everystep = false,
    saveat = 0.1, abstol = 1.0e-5, reltol = 1.0e-3
)

println(sol.u[end])

#[0.011503,0.939809,0.00312214,0.00155873,0.0172325,0.0577331,0.237757,0.00134921,0.000238022,4.19916e-5,7.40824e-6,1.307e-6,0.0621091,1.24463,0.0483949,199.901,137.457,0.0177237,0.132583]

Random.seed!(100)
@time integrator = init(
    quick_prob, SRIW1(), dt = (1 / 2)^(18), progress_steps = Int(1.0e5),
    adaptivealg = :RSwM3, progress = false, qmax = 4, save_everystep = false,
    saveat = 0.1, abstol = 1.0e-5, reltol = 1.0e-3
)

@time StochasticDiffEq.solve!(integrator)

@code_warntype StochasticDiffEq.generate_tildes(integrator, 1.0, 1.0, 1.0)
@code_warntype StochasticDiffEq.perform_step!(integrator)
@code_llvm DiffEqBase.isinplace(integrator.noise)
function test_pop(integrator)
    return b2 = integrator.opts.beta2
end
@code_warntype test_pop(integrator)
