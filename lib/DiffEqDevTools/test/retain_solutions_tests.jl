using DiffEqDevTools, StochasticDiffEq, Random, Test

# The Monte-Carlo `test_convergence` summarises by default so that large studies do not
# hold every trajectory of every step size at once. That must change only what is kept,
# never what is measured, so compare the statistics against `retain_solutions = true`
# exactly. It also pins the default, which is the breaking half of the change.

linear_analytic(u0, p, t, W) = @.(u0 * exp(0.63155t + 0.87W))
f_linear(du, u, p, t) = @.(du = 1.01 * u)
σ_linear(du, u, p, t) = @.(du = 0.87 * u)
prob = SDEProblem(
    SDEFunction(f_linear, σ_linear, analytic = linear_analytic), [1 / 2], (0.0, 1.0)
)
dts = 1 .// 2 .^ (7:-1:4)
trajectories = 400

function run(retain)
    Random.seed!(100)
    return test_convergence(
        dts, prob, EM(), EnsembleSerial(); save_everystep = false,
        trajectories = trajectories, weak_timeseries_errors = false,
        retain_solutions = retain
    )
end

kept = run(true)
dropped = run(false)
defaulted = let
    Random.seed!(100)
    test_convergence(
        dts, prob, EM(), EnsembleSerial(); save_everystep = false,
        trajectories = trajectories, weak_timeseries_errors = false
    )
end

@testset "retain_solutions preserves the statistics" begin
    @test keys(kept.𝒪est) == keys(dropped.𝒪est)
    for k in keys(kept.𝒪est)
        @test kept.𝒪est[k] == dropped.𝒪est[k]
    end
    @test keys(kept.errors) == keys(dropped.errors)
    for k in keys(kept.errors)
        @test kept.errors[k] == dropped.errors[k]
    end
    for i in eachindex(dts)
        @test kept.solutions[i].weak_errors == dropped.solutions[i].weak_errors
        @test kept.solutions[i].error_means == dropped.solutions[i].error_means
        @test kept.solutions[i].errors == dropped.solutions[i].errors
    end
end

@testset "summarising is the default" begin
    for i in eachindex(dts)
        @test length(defaulted.solutions[i].u) == 1
        @test defaulted.solutions[i].errors == kept.solutions[i].errors
    end
    for k in keys(kept.𝒪est)
        @test defaulted.𝒪est[k] == kept.𝒪est[k]
    end
end

@testset "retain_solutions drops the trajectories" begin
    for i in eachindex(dts)
        @test length(kept.solutions[i].u) == trajectories
        @test length(dropped.solutions[i].u) == 1
    end
    @test Base.summarysize(dropped) < Base.summarysize(kept) / 10
end

@testset "retain_solutions is ignored when expected_value is given" begin
    # The expected_value path takes a reducing output_func and never builds an
    # EnsembleTestSolution, so there is nothing to drop; it must still run and agree.
    ev = 0.5 * exp(1.01)
    ensemble_prob = EnsembleProblem(prob; output_func = (sol, ctx) -> (sol.u[end][1], false))
    Random.seed!(100)
    a = test_convergence(
        dts, ensemble_prob, EM(), EnsembleSerial(); save_everystep = false,
        trajectories = trajectories, expected_value = ev, retain_solutions = false
    )
    Random.seed!(100)
    b = test_convergence(
        dts, ensemble_prob, EM(), EnsembleSerial(); save_everystep = false,
        trajectories = trajectories, expected_value = ev
    )
    @test a.𝒪est[:weak_final] == b.𝒪est[:weak_final]
    @test length(a.solutions[1].u) == trajectories
end
