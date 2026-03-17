using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqROCK
        using Test

        @test SROCK1() isa StochasticDiffEqAlgorithm
        @test SROCK2() isa StochasticDiffEqAlgorithm
        @test KomBurSROCK2() isa StochasticDiffEqAlgorithm
        @test SROCKC2() isa StochasticDiffEqAlgorithm
        @test SROCKEM() isa StochasticDiffEqAlgorithm
        @test SKSROCK() isa StochasticDiffEqAlgorithm
        @test TangXiaoSROCK2() isa StochasticDiffEqAlgorithm
    end
end
