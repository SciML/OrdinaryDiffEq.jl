using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqIIF
        using Test

        @test IIF1M() isa IIF1M
        @test IIF2M() isa IIF2M
        @test IIF1Mil() isa IIF1Mil
        @test StochasticDiffEqIIF.alg_order(IIF1M()) == 1 // 2
        @test StochasticDiffEqIIF.alg_order(IIF2M()) == 1 // 2
        @test StochasticDiffEqIIF.alg_order(IIF1Mil()) == 1 // 1
    end
end
