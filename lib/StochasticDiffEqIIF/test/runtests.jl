using StochasticDiffEqIIF
using Test

@testset "StochasticDiffEqIIF" begin
    @testset "Module loading" begin
        @test true
    end

    @testset "Algorithm constructors" begin
        @test IIF1M() isa StochasticDiffEqIIF.IIF1M
        @test IIF2M() isa StochasticDiffEqIIF.IIF2M
        @test IIF1Mil() isa StochasticDiffEqIIF.IIF1Mil
    end

    @testset "Algorithm orders" begin
        @test StochasticDiffEqIIF.alg_order(IIF1M()) == 1 // 2
        @test StochasticDiffEqIIF.alg_order(IIF2M()) == 1 // 2
        @test StochasticDiffEqIIF.alg_order(IIF1Mil()) == 1 // 1
    end
end
