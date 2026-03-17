using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqHighOrder
        using Test

        @test SRI() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRIW1() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRIW2() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRI() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRI2() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA1() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA2() isa StochasticDiffEqAdaptiveAlgorithm
        @test SRA3() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRA() isa StochasticDiffEqAdaptiveAlgorithm
        @test SOSRA2() isa StochasticDiffEqAdaptiveAlgorithm
    end
end
