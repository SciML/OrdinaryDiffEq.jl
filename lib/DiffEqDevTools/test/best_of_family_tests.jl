using OrdinaryDiffEq, DiffEqDevTools, Test
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear

@testset "Best-of-Family" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => Tsit5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => BS3(), :tags => [:rk, :third_order]),
        Dict{Symbol, Any}(:alg => Euler(), :tags => [:euler, :first_order]),
    ]

    wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2)

    @testset "wp_area" begin
        area = wp_area(wp_set[1])
        @test isfinite(area)
    end

    @testset "best_by_tag" begin
        best_rk = best_by_tag(wp_set, :rk; n = 2)
        @test length(best_rk) == 2

        best_euler = best_by_tag(wp_set, :euler; n = 1)
        @test length(best_euler) == 1
        @test best_euler.names == ["Euler"]

        # No match
        empty = best_by_tag(wp_set, :nonexistent)
        @test length(empty) == 0
    end

    @testset "best_of_families" begin
        best = best_of_families(wp_set, [:rk, :euler]; n = 1)
        @test length(best) == 2
        @test "Euler" in best.names

        # Error on all-empty
        @test_throws ArgumentError best_of_families(wp_set, [:nonexistent])
    end

    @testset "unknown metric" begin
        @test_throws ArgumentError best_by_tag(wp_set, :rk; metric = :unknown)
    end
end
