using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqMilstein
        using Test

        @test RKMilGeneral() isa StochasticDiffEqAdaptiveAlgorithm
        @test WangLi3SMil_A() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_B() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_C() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_D() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_E() isa StochasticDiffEqAlgorithm
        @test WangLi3SMil_F() isa StochasticDiffEqAlgorithm
    end
end
