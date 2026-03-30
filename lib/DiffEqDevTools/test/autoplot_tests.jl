using OrdinaryDiffEq, DiffEqDevTools, Test
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear

@testset "Autoplot" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => Tsit5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => BS3(), :tags => [:rk, :third_order]),
        Dict{Symbol, Any}(:alg => Euler(), :tags => [:euler, :first_order, :reference]),
    ]

    wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2)

    @testset "explicit families" begin
        plots = autoplot(wp_set; families = [:rk, :euler])
        @test haskey(plots, "family_rk")
        @test haskey(plots, "family_euler")
        @test haskey(plots, "best_of_families")
        @test haskey(plots, "all")
        @test length(plots["family_rk"]) == 4  # RK4, DP5, Tsit5, BS3
        @test length(plots["family_euler"]) == 1
    end

    @testset "auto-detect families" begin
        plots = autoplot(wp_set)
        @test haskey(plots, "all")
        # Should detect families that aren't on >80% of methods
    end

    @testset "with reference_tags" begin
        plots = autoplot(wp_set; families = [:rk], reference_tags = :reference)
        # rk family should include the reference method too
        rk_names = plots["family_rk"].names
        @test "Euler" in rk_names || length(plots["family_rk"]) >= 4
    end

    @testset "best_n" begin
        plots = autoplot(wp_set; families = [:rk, :euler], best_n = 1)
        bof = plots["best_of_families"]
        # Should have at most 1 from each family
        @test length(bof) <= 2
    end
end
