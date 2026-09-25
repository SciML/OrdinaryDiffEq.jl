using OrdinaryDiffEq, DiffEqDevTools, Test
using OrdinaryDiffEqLowOrderRK: RK4
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear
abstols = 1 ./ 10 .^ (3:6)
reltols = 1 ./ 10 .^ (3:6)

@testset "no timeout leaves every point measured" begin
    wp = WorkPrecision(
        prob, RK4(), abstols, reltols; name = "RK4", dt = 1 / 2^4, numruns = 2,
        timeout = nothing
    )
    @test !any(isnan, wp.times)
end

@testset "an exceeded budget records NaN" begin
    wp = @test_logs (:warn,) (:warn,) (:warn,) (:warn,) WorkPrecision(
        prob, RK4(), abstols, reltols; name = "RK4", dt = 1 / 2^4, numruns = 2,
        timeout = 0.0
    )
    @test all(isnan, wp.times)
    @test all(isnan, wp.errors.final)
    # the solve that blew the budget still ran, so its statistics are kept
    @test all(!isnothing, wp.stats)
end

@testset "WorkPrecisionSet forwards the budget" begin
    setups = [Dict{Symbol, Any}(:alg => RK4())]
    wp_set = @test_logs (:warn,) (:warn,) (:warn,) (:warn,) WorkPrecisionSet(
        prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2, timeout = 0.0
    )
    @test all(isnan, wp_set[1].times)
end
