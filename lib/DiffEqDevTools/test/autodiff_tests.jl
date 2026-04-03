using DiffEqDevTools, Test

# Mock AD backend types for testing
struct AutoForwardDiff end
struct AutoFiniteDiff end
struct ManualJac end

@testset "AutoDiff Helpers" begin
    @testset "ad_backend_name" begin
        @test DiffEqDevTools.ad_backend_name(AutoForwardDiff()) == "ForwardDiff"
        @test DiffEqDevTools.ad_backend_name(AutoFiniteDiff()) == "FiniteDiff"
        @test DiffEqDevTools.ad_backend_name(ManualJac()) == "ManualJac"
    end

    @testset "with_autodiff_variants count" begin
        setups = [
            Dict{Symbol, Any}(:alg => nothing, :tags => [:rk]),
            Dict{Symbol, Any}(:alg => nothing),
        ]
        expanded = with_autodiff_variants(setups;
            ad_backends = [AutoForwardDiff(), AutoFiniteDiff()])
        # 2 originals + 2x2 variants = 6
        @test length(expanded) == 6
    end

    @testset "original gets default tag" begin
        setups = [Dict{Symbol, Any}(:alg => nothing, :tags => [:rk])]
        expanded = with_autodiff_variants(setups; ad_backends = [AutoForwardDiff()])
        @test :autodiff_default in expanded[1][:tags]
        @test :rk in expanded[1][:tags]
    end

    @testset "variant gets backend tag and autodiff key" begin
        setups = [Dict{Symbol, Any}(:alg => nothing, :tags => [:rk])]
        expanded = with_autodiff_variants(setups; ad_backends = [AutoForwardDiff()])
        @test expanded[2][:autodiff] isa AutoForwardDiff
        @test :rk in expanded[2][:tags]
    end

    @testset "custom tag_prefix" begin
        setups = [Dict{Symbol, Any}(:alg => nothing)]
        expanded = with_autodiff_variants(setups;
            ad_backends = [AutoForwardDiff()], tag_prefix = :ad)
        @test :ad_default in expanded[1][:tags]
    end

    @testset "no mutation of original" begin
        original_setup = Dict{Symbol, Any}(:alg => nothing, :tags => [:rk])
        setups = [original_setup]
        expanded = with_autodiff_variants(setups; ad_backends = [AutoForwardDiff()])
        @test :autodiff_default ∉ original_setup[:tags]
        @test length(original_setup[:tags]) == 1
    end
end
