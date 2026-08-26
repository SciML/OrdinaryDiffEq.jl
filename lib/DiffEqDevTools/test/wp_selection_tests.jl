using OrdinaryDiffEq, DiffEqDevTools, Test
using OrdinaryDiffEqLowOrderRK: Euler, BS3, DP5, RK4
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqRosenbrock: Rosenbrock23
using OrdinaryDiffEqSDIRK: TRBDF2
using DiffEqDevTools: _default_name
using ODEProblemLibrary: prob_ode_linear

prob = prob_ode_linear
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

@testset "wp_area" begin
    @test all(isfinite, wp_area.(wp_set.wps))

    # A curve with fewer than two usable points has no area to integrate.
    lonely = deepcopy(wp_set[1])
    lonely.times = [NaN; fill(NaN, length(lonely.times) - 2); 1.0]
    @test wp_area(lonely) == Inf
end

@testset "best_by_tag" begin
    @test length(best_by_tag(wp_set, :rk; n = 2)) == 2
    @test best_by_tag(wp_set, :euler).names == ["Euler"]
    @test length(@test_logs (:warn,) best_by_tag(wp_set, :nonexistent)) == 0
    @test_throws ArgumentError best_by_tag(wp_set, :rk; metric = :unknown)

    # `n` truncates one ranking rather than re-ranking each time
    @test best_by_tag(wp_set, :rk; n = 2).names == best_by_tag(wp_set, :rk; n = 4).names[1:2]

    # Ranking prefers the method that produced usable results at every tolerance.
    crippled = deepcopy(wp_set)
    crippled.wps[2].times = [crippled.wps[2].times[1]; fill(NaN, length(abstols) - 1)]
    @test "DP5" ∉ best_by_tag(crippled, :rk; n = 3).names
end

@testset "best_of_families" begin
    best = best_of_families(wp_set, [:rk, :euler]; n = 1)
    @test length(best) == 2
    @test "Euler" in best.names
    @test_throws ArgumentError best_of_families(wp_set, [:nonexistent])

    # A method in two families is included once
    overlapping = best_of_families(wp_set, [:rk, :fifth_order]; n = 5)
    @test length(overlapping) == 4
    @test allunique(overlapping.names)
end

@testset "autoplot" begin
    plots = autoplot(wp_set; families = [:rk, :euler])
    @test sort(collect(keys(plots))) == ["all", "best_of_families", "family_euler", "family_rk"]
    @test length(plots["family_rk"]) == 4
    @test length(plots["family_euler"]) == 1
    @test plots["all"] === wp_set
    @test length(plots["best_of_families"]) == 3  # best 2 of :rk, and the lone :euler

    @testset "reference methods join every family" begin
        plots = autoplot(wp_set; families = [:rk], reference_tags = :reference)
        @test plots["family_rk"].names == ["RK4", "DP5", "Tsit5", "BS3", "Euler"]
        @test "Euler" in plots["best_of_families"].names
    end

    @testset "auto-detected families" begin
        # :rk is on 4/5 entries, :fifth_order on 2/5; :explicit is on every entry
        # and is therefore dropped.
        @test DiffEqDevTools._auto_detect_families(wp_set) ==
            filter(!=(:explicit), unique_tags(wp_set))
        uniform = filter_by_tags(wp_set, :rk)
        @test :rk ∉ DiffEqDevTools._auto_detect_families(uniform)
        @test haskey(autoplot(wp_set), "family_fifth_order")
    end

    @testset "families with no entries are skipped" begin
        plots = autoplot(wp_set; families = [:rk, :nonexistent])
        @test !haskey(plots, "family_nonexistent")
    end
end

@testset "with_autodiff_variants" begin
    original = Dict{Symbol, Any}(:alg => Rosenbrock23(), :tags => [:rosenbrock])
    expanded = with_autodiff_variants(
        [original, Dict{Symbol, Any}(:alg => TRBDF2()), Dict{Symbol, Any}(:alg => Tsit5())];
        ad_backends = [AutoForwardDiff(), AutoFiniteDiff()]
    )

    @test DiffEqDevTools.ad_backend_name(AutoForwardDiff()) == "ForwardDiff"

    # two implicit setups expand to three entries each, Tsit5 takes no autodiff
    @test length(expanded) == 7
    @test expanded[1][:tags] == [:rosenbrock, :autodiff_default]
    @test expanded[2][:tags] == [:rosenbrock, :autodiff_forwarddiff]
    @test expanded[2][:alg].autodiff isa AutoForwardDiff
    @test expanded[3][:alg].autodiff isa AutoFiniteDiff
    @test expanded[3][:name] == "Rosenbrock23 (FiniteDiff)"
    @test expanded[7][:tags] == [:autodiff_default]
    @test _default_name(expanded[7][:alg]) == "Tsit5"
    @test !haskey(expanded[1], :name)

    @test original[:tags] == [:rosenbrock]  # inputs are not mutated

    @test with_autodiff_variants(
        [original]; ad_backends = [AutoForwardDiff()], tag_prefix = :ad
    )[1][:tags] == [:rosenbrock, :ad_default]

    @testset "the variants reach the benchmark" begin
        ad_set = WorkPrecisionSet(
            prob, abstols, reltols,
            with_autodiff_variants(
                [Dict{Symbol, Any}(:alg => Rosenbrock23())];
                ad_backends = [AutoFiniteDiff()]
            ); numruns = 2
        )
        @test ad_set.names == ["Rosenbrock23", "Rosenbrock23 (FiniteDiff)"]
        @test get_tags(ad_set) == [
            [:order_2, :adaptive, :implicit, :rosenbrock, :autodiff_default],
            [:order_2, :adaptive, :implicit, :rosenbrock, :autodiff_finitediff],
        ]
        @test ad_set.setups[2][:alg].autodiff isa AutoFiniteDiff
    end
end
