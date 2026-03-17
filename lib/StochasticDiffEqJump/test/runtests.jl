using StochasticDiffEqJump
using Test

@testset "StochasticDiffEqJump" begin
    @testset "Module loading" begin
        @test @isdefined(TauLeaping)
        @test @isdefined(CaoTauLeaping)
        @test @isdefined(ImplicitTauLeaping)
        @test @isdefined(ThetaTrapezoidalTauLeaping)
    end

    @testset "Algorithm construction" begin
        @test TauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test CaoTauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test ImplicitTauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test ThetaTrapezoidalTauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test ThetaTrapezoidalTauLeaping(theta = 0.7).theta == 0.7
    end
end
