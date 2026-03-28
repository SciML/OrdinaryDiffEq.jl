using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear

@testset "Tags in WorkPrecision" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    tgs = Vector{Symbol}([:rk, :fifth_order])
    wp = WorkPrecision(
        prob, DP5(), abstols, reltols;
        name = "DP5", tags = tgs
    )
    @test wp.tags == [:rk, :fifth_order]

    wp_notags = WorkPrecision(prob, DP5(), abstols, reltols; name = "DP5")
    @test isempty(wp_notags.tags)
end

@testset "Tags in WorkPrecisionSet via setups" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => BS3(), :tags => [:rk, :third_order]),
        Dict{Symbol, Any}(:alg => Euler(), :tags => [:rk, :first_order, :reference]),
        Dict{Symbol, Any}(:alg => Midpoint()),
    ]

    wp_set = WorkPrecisionSet(
        prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2
    )

    @test wp_set[1].tags == [:rk, :fourth_order]
    @test wp_set[2].tags == [:rk, :fifth_order]
    @test wp_set[3].tags == [:rk, :third_order]
    @test wp_set[4].tags == [:rk, :first_order, :reference]
    @test isempty(wp_set[5].tags)
end

@testset "filter_by_tags" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => BS3(), :tags => [:rk, :third_order]),
        Dict{Symbol, Any}(:alg => Euler(), :tags => [:euler, :first_order, :reference]),
        Dict{Symbol, Any}(:alg => Midpoint()),
    ]

    wp_set = WorkPrecisionSet(
        prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2
    )

    # Filter by single tag
    rk_only = filter_by_tags(wp_set, :rk)
    @test length(rk_only) == 3
    @test Set(rk_only.names) == Set(["RK4", "DP5", "BS3"])

    # Filter by multiple tags (AND logic)
    rk_5th = filter_by_tags(wp_set, :rk, :fifth_order)
    @test length(rk_5th) == 1
    @test rk_5th.names == ["DP5"]

    # Filter by tag that matches nothing
    empty_result = filter_by_tags(wp_set, :nonexistent)
    @test length(empty_result) == 0

    # Empty tags returns original
    same = filter_by_tags(wp_set)
    @test length(same) == length(wp_set)
end

@testset "exclude_by_tags" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => Euler(), :tags => [:euler, :reference]),
        Dict{Symbol, Any}(:alg => Midpoint()),
    ]

    wp_set = WorkPrecisionSet(
        prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2
    )

    # Exclude reference
    no_ref = exclude_by_tags(wp_set, :reference)
    @test length(no_ref) == 3
    @test "Euler" ∉ no_ref.names

    # Exclude multiple (OR logic)
    no_rk_or_ref = exclude_by_tags(wp_set, :rk, :reference)
    @test length(no_rk_or_ref) == 1
    @test no_rk_or_ref.names == ["Midpoint"]

    # Empty exclude returns original
    same = exclude_by_tags(wp_set)
    @test length(same) == length(wp_set)
end

@testset "unique_tags and get_tags" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk, :fourth_order]),
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk, :fifth_order]),
        Dict{Symbol, Any}(:alg => Euler(), :tags => [:euler, :reference]),
    ]

    wp_set = WorkPrecisionSet(
        prob, abstols, reltols, setups; dt = 1 / 2^4, numruns = 2
    )

    tags = get_tags(wp_set)
    @test tags[1] == [:rk, :fourth_order]
    @test tags[2] == [:rk, :fifth_order]
    @test tags[3] == [:euler, :reference]

    utags = unique_tags(wp_set)
    @test Set(utags) == Set([:rk, :fourth_order, :fifth_order, :euler, :reference])
end

@testset "merge_wp_sets" begin
    abstols = 1 ./ 10 .^ (3:6)
    reltols = 1 ./ 10 .^ (3:6)

    setups1 = [
        Dict{Symbol, Any}(:alg => RK4(), :tags => [:rk]),
    ]
    setups2 = [
        Dict{Symbol, Any}(:alg => DP5(), :tags => [:rk]),
        Dict{Symbol, Any}(:alg => BS3(), :tags => [:rk]),
    ]

    wp1 = WorkPrecisionSet(prob, abstols, reltols, setups1; dt = 1 / 2^4, numruns = 2)
    wp2 = WorkPrecisionSet(prob, abstols, reltols, setups2; dt = 1 / 2^4, numruns = 2)

    merged = merge_wp_sets(wp1, wp2)
    @test length(merged) == 3
    @test merged.names == ["RK4", "DP5", "BS3"]
end
