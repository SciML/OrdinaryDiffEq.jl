using StochasticDiffEqImplicit
using Test

@testset "StochasticDiffEqImplicit" begin
    @testset "Module loads" begin
        @test true
    end

    @testset "Algorithm construction" begin
        @test ImplicitEM() isa StochasticDiffEqImplicit.StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitEulerHeun() isa StochasticDiffEqImplicit.StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitRKMil() isa StochasticDiffEqImplicit.StochasticDiffEqNewtonAdaptiveAlgorithm
        @test STrapezoid() isa StochasticDiffEqImplicit.ImplicitEM
        @test SImplicitMidpoint() isa StochasticDiffEqImplicit.ImplicitEM
        @test ISSEM() isa StochasticDiffEqImplicit.StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ISSEulerHeun() isa StochasticDiffEqImplicit.StochasticDiffEqNewtonAdaptiveAlgorithm
        @test SKenCarp() isa StochasticDiffEqImplicit.StochasticDiffEqNewtonAdaptiveAlgorithm
    end
end
