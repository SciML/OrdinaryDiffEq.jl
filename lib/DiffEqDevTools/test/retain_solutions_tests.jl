using StochasticDiffEq, DiffEqDevTools, Random, Test

# The Monte-Carlo method of `test_convergence` used to keep one full solution object per
# trajectory per step size, which is orders of magnitude more than the error estimates
# need and is what put the weak-convergence groups past a 16 GB runner. Trajectories are
# now reduced by default: an `output_func` turns each into a `ConvergenceTrajectory` as
# it is solved, and the ensemble is stripped to one representative once the errors have
# been computed. The statistics must be identical either way.

f_iip(du, u, p, t) = @.(du = 1.01 * u)
σ_iip(du, u, p, t) = @.(du = 0.87 * u)
linear_analytic(u0, p, t, W) = @.(u0 * exp(0.63155t + 0.87W))
prob = SDEProblem(
    SDEFunction(f_iip, σ_iip, analytic = linear_analytic), [1 / 2], (0.0, 1.0)
)

dts = 1 .// 2 .^ (5:-1:3)
ntraj = 400

function run_sim(; kwargs...)
    Random.seed!(100)
    return test_convergence(
        dts, prob, EM(); save_everystep = false, trajectories = ntraj,
        weak_timeseries_errors = false, kwargs...
    )
end

reduced = run_sim()                            # the new default
full = run_sim(retain_solutions = true)

@testset "statistics are unchanged" begin
    @test keys(full.𝒪est) == keys(reduced.𝒪est)
    for k in keys(full.𝒪est)
        @test full.𝒪est[k] ≈ reduced.𝒪est[k]
    end
    @test keys(full.errors) == keys(reduced.errors)
    for k in keys(full.errors)
        @test collect(full.errors[k]) ≈ collect(reduced.errors[k])
    end
    for i in eachindex(dts)
        @test full.solutions[i].weak_errors == reduced.solutions[i].weak_errors
        @test full.solutions[i].error_means == reduced.solutions[i].error_means
        # the per-trajectory error vectors are computed from the full ensemble either way
        for k in keys(full.solutions[i].errors)
            @test full.solutions[i].errors[k] ≈ reduced.solutions[i].errors[k]
            @test length(reduced.solutions[i].errors[k]) == ntraj
        end
    end
end

@testset "trajectories are dropped by default" begin
    for i in eachindex(dts)
        @test length(reduced.solutions[i].u) == 1
        @test length(full.solutions[i].u) == ntraj
    end
    @test reduced.solutions[1].u[1] isa ConvergenceTrajectory
    @test full.solutions[1].u[1] isa SciMLBase.AbstractRODESolution
    # the representative keeps the endpoints the weak error is formed from
    @test reduced.solutions[1].u[1].u[end] == full.solutions[1].u[1].u[end]
    @test reduced.solutions[1].u[1].u_analytic[end] == full.solutions[1].u[1].u_analytic[end]
    @test Base.summarysize(reduced.solutions) < Base.summarysize(full.solutions) / 20
end

@testset "reduction is skipped where it cannot apply" begin
    # these need the full timeseries, so the solutions are retained rather than the
    # call failing
    for kwargs in ((; weak_timeseries_errors = true), (; weak_dense_errors = true))
        sim = run_sim(; kwargs...)
        @test length(sim.solutions[1].u) == ntraj
    end
    # a caller-supplied EnsembleProblem owns its own output_func
    Random.seed!(100)
    sim = test_convergence(
        dts, EnsembleProblem(prob), EM(); save_everystep = false,
        trajectories = ntraj, weak_timeseries_errors = false
    )
    @test length(sim.solutions[1].u) == ntraj

    # `expected_value` averages the trajectory values themselves, so it is used with an
    # EnsembleProblem that reduces to the quantity being averaged (as PL1WM.jl does).
    # That path must keep every reduced value, not one representative.
    Random.seed!(100)
    ens = EnsembleProblem(prob; output_func = (sol, ctx) -> (sol.u[end][1], false))
    sim = test_convergence(
        dts, ens, EM(); save_everystep = false, save_start = false,
        trajectories = ntraj, weak_timeseries_errors = false,
        expected_value = 1 / 2 * exp(1.01)
    )
    @test length(sim.solutions[1].u) == ntraj
    @test sim.solutions[1].u[1] isa Float64
    @test haskey(sim.𝒪est, :weak_final)
end
