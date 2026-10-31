using OrdinaryDiffEq, DiffEqDevTools, Test
using OrdinaryDiffEqLowOrderRK: DP5, RK4
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear
abstols = 1 ./ 10 .^ (3:6)
reltols = 1 ./ 10 .^ (3:6)
setups = [Dict{Symbol, Any}(:alg => DP5())]

@testset "one error estimate (the default)" begin
    wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; numruns = 2)
    @test available_errors(wp_set) == [:final]
    @test propertynames(wp_set[1].errors) == (:final,)
end

@testset "several error estimates from one run" begin
    wp_set = WorkPrecisionSet(
        prob, abstols, reltols, setups; numruns = 2,
        error_estimates = [:final, :l2, :L2]
    )
    @test available_errors(wp_set) == [:final, :l2, :L2]
    for e in available_errors(wp_set)
        @test all(isfinite, getproperty(wp_set[1].errors, e))
    end
end

@testset "positional constructors without the new fields" begin
    wp = WorkPrecision(prob, DP5(), abstols, reltols; name = "DP5", numruns = 2)
    bare = WorkPrecision(
        prob, abstols, reltols, wp.errors, wp.times, wp.dts, wp.stats, "DP5", :final,
        length(abstols)
    )
    @test isempty(bare.tags)

    wp_set = WorkPrecisionSet(
        [wp], 1, abstols, reltols, prob, setups, ["DP5"], :final, nothing
    )
    @test available_errors(wp_set) == [:final]
end

@testset "tolerances that fail keep the successful ones plottable" begin
    # Regression test: a failed solve records every error key while a successful one
    # records only the requested estimate, which StructArrays 0.7 rejects unless the
    # key sets are reconciled.
    wp = WorkPrecision(
        prob, RK4(), [1.0e-3, 1.0e-14], [1.0e-3, 1.0e-14];
        name = "RK4", numruns = 2, maxiters = 300
    )
    @test !isnan(wp.times[1])
    @test isnan(wp.times[2])
    @test isfinite(wp.errors.final[1])
    @test isnan(wp.errors.final[2])
end
