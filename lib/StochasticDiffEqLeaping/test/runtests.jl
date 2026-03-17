using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqLeaping
        using Test

        @test TauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test CaoTauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test ImplicitTauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test ThetaTrapezoidalTauLeaping() isa StochasticDiffEqJumpAdaptiveAlgorithm
        @test ThetaTrapezoidalTauLeaping(theta = 0.7).theta == 0.7
    end
end
