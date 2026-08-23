using OrdinaryDiffEq, DiffEqDevTools, Test
using OrdinaryDiffEqLowOrderRK: Euler, Midpoint, BS3, DP5, RK4
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear
abstols = 1 ./ 10 .^ (3:6)
reltols = 1 ./ 10 .^ (3:6)

setups = [
    Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
    Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
    Dict{Symbol, Any}(:alg => BS3(), :tags => [:rk, :third_order]),
    Dict{Symbol, Any}(:alg => Euler(), :tags => [:euler, :first_order, :reference]),
    Dict{Symbol, Any}(:alg => Midpoint()),
]
wp_set = WorkPrecisionSet(prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2)

@testset "tags reach the WorkPrecision entries" begin
    @test get_tags(wp_set) == [
        [:rk, :fourth_order], [:rk, :fifth_order], [:rk, :third_order],
        [:euler, :first_order, :reference], Symbol[],
    ]
    @test unique_tags(wp_set) ==
        sort([:rk, :fourth_order, :fifth_order, :third_order, :euler, :first_order, :reference])

    wp = WorkPrecision(
        prob, DP5(), abstols, reltols; name = "DP5", numruns = 2,
        tags = [:rk, :fifth_order]
    )
    @test wp.tags == [:rk, :fifth_order]
    @test isempty(WorkPrecision(prob, DP5(), abstols, reltols; numruns = 2).tags)
end

@testset "filter_by_tags" begin
    for (tags, expected) in [
            ((:rk,), ["RK4", "DP5", "BS3"]),
            ((:rk, :fifth_order), ["DP5"]),
            ((:reference,), ["Euler"]),
        ]
        subset = filter_by_tags(wp_set, tags...)
        @test subset.names == expected
        @test [wp.name for wp in subset.wps] == expected
        # setups stay aligned with the entries that survived the filter
        @test [string(nameof(typeof(s[:alg]))) for s in subset.setups] == expected
    end

    @test length(@test_logs (:warn,) filter_by_tags(wp_set, :nonexistent)) == 0
    @test length(filter_by_tags(wp_set)) == length(wp_set)
end

@testset "exclude_by_tags" begin
    @test exclude_by_tags(wp_set, :reference).names == ["RK4", "DP5", "BS3", "Midpoint"]
    @test exclude_by_tags(wp_set, :rk, :reference).names == ["Midpoint"]
    @test length(exclude_by_tags(wp_set)) == length(wp_set)
end

@testset "merge_wp_sets" begin
    merged = merge_wp_sets(filter_by_tags(wp_set, :euler), filter_by_tags(wp_set, :rk))
    @test length(merged) == 4
    @test merged.names == ["Euler", "RK4", "DP5", "BS3"]
    @test get_tags(merged)[1] == [:euler, :first_order, :reference]
    @test available_errors(merged) == available_errors(wp_set)
    @test_throws ArgumentError merge_wp_sets()
end
