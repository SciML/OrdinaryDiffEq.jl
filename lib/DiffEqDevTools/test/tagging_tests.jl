using OrdinaryDiffEq, DiffEqDevTools, Test
using OrdinaryDiffEqFIRK: RadauIIA5
using OrdinaryDiffEqLowOrderRK: Euler, Midpoint, BS3, DP5, RK4
using OrdinaryDiffEqRosenbrock: Rosenbrock23
using OrdinaryDiffEqSDIRK: KenCarp4, TRBDF2
using StochasticDiffEqHighOrder: SRIW1
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

@testset "automatic algorithm tags" begin
    for (alg, expected) in [
            (Euler(), [:order_1, :fixed_step, :explicit]),
            (DP5(), [:order_5, :adaptive, :explicit]),
            (
                Rosenbrock23(),
                [:order_2, :adaptive, :implicit, :rosenbrock],
            ),
            (
                TRBDF2(),
                [:order_2, :adaptive, :implicit, :sdirk, :esdirk],
            ),
            (
                KenCarp4(),
                [:order_4, :adaptive, :implicit, :sdirk, :esdirk, :split],
            ),
            (RadauIIA5(), [:order_5, :adaptive, :implicit, :firk]),
            (SRIW1(), [:order_3_2, :adaptive, :sde]),
        ]
        @test auto_tags(alg) == expected
    end
end

DiffEqDevTools.tag_kind(::Val{:package_family}) = :family

@testset "tag semantics" begin
    for (tag, expected) in [
            (:rk, :family),
            (:rosenbrock, :family),
            (:order_5, :trait),
            (:adaptive, :trait),
            (:reference, :role),
            (:sundials, :provider),
            (:dense_output, :variant),
            (:stiff, :domain),
            (:package_family, :family),
            (:benchmark_specific, :unknown),
        ]
        @test tag_kind(tag) == expected
    end
end

@testset "tags reach the WorkPrecision entries" begin
    @test get_tags(wp_set) == [
        [:order_4, :adaptive, :explicit, :rk, :fourth_order],
        [:order_5, :adaptive, :explicit, :rk, :fifth_order],
        [:order_3, :adaptive, :explicit, :rk, :third_order],
        [:order_1, :fixed_step, :explicit, :euler, :first_order, :reference],
        [:order_2, :adaptive, :explicit],
    ]
    @test unique_tags(wp_set) ==
        sort(
        [
            :adaptive, :euler, :explicit, :fifth_order, :first_order, :fixed_step,
            :fourth_order, :order_1, :order_2, :order_3, :order_4, :order_5, :reference,
            :rk, :third_order,
        ]
    )

    wp = WorkPrecision(
        prob, DP5(), abstols, reltols; name = "DP5", numruns = 2,
        tags = [:rk, :fifth_order]
    )
    @test wp.tags == [:order_5, :adaptive, :explicit, :rk, :fifth_order]
    @test WorkPrecision(prob, DP5(), abstols, reltols; numruns = 2).tags ==
        [:order_5, :adaptive, :explicit]
    @test WorkPrecision(
        prob, DP5(), abstols, reltols; numruns = 2,
        tags = [:manual], auto_tags = false
    ).tags == [:manual]
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
    @test get_tags(merged)[1] ==
        [:order_1, :fixed_step, :explicit, :euler, :first_order, :reference]
    @test available_errors(merged) == available_errors(wp_set)
    @test_throws ArgumentError merge_wp_sets()
end
