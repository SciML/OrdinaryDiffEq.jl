using OrdinaryDiffEq, DiffEqDevTools, Test
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear

@testset "Multi-error-mode" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    @testset "Backward compatible 9-arg constructor" begin
        wp = WorkPrecision(prob, DP5(), abstols, reltols; name = "DP5", numruns = 2)
        wps = WorkPrecisionSet([wp], 1, abstols, reltols, prob, [], ["DP5"], :final, nothing)
        @test wps.active_error_estimates == [:final]
        @test available_errors(wps) == [:final]
    end

    @testset "Single error_estimate (default)" begin
        setups = [Dict{Symbol, Any}(:alg => DP5())]
        wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; numruns = 2)
        @test available_errors(wp_set) == [:final]
    end

    @testset "Multi error_estimates" begin
        setups = [Dict{Symbol, Any}(:alg => DP5())]
        wp_set = WorkPrecisionSet(prob, abstols, reltols, setups;
            numruns = 2, error_estimates = [:final, :l2])
        @test Set(available_errors(wp_set)) == Set([:final, :l2])
        # l2 errors should be computed (timeseries_errors was enabled)
        @test :l2 in propertynames(wp_set[1].errors)
    end
end
