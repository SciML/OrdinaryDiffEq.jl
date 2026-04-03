using OrdinaryDiffEq, DiffEqDevTools, Test
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear

@testset "Timeout Tests" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    @testset "timeout=nothing preserves behavior" begin
        wp = WorkPrecision(prob, RK4(), abstols, reltols;
            name = "RK4", dt = 1 / 2^4, numruns = 2, timeout = nothing)
        @test !any(isnan, wp.times)
    end

    @testset "timeout=0.0 causes NaN" begin
        wp = @test_warn r"Timeout" WorkPrecision(prob, RK4(), abstols, reltols;
            name = "RK4", dt = 1 / 2^4, numruns = 2, timeout = 0.0)
        @test all(isnan, wp.times)
    end

    @testset "WorkPrecisionSet passes timeout" begin
        setups = [Dict{Symbol, Any}(:alg => RK4())]
        wp_set = @test_warn r"Timeout" WorkPrecisionSet(prob, abstols, reltols, setups;
            dt = 1 / 2^4, numruns = 2, timeout = 0.0)
        @test all(isnan, wp_set[1].times)
    end
end
