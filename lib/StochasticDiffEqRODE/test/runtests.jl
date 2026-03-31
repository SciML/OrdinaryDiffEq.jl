using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqRODE
        using Test

        @test RandomEM() isa StochasticDiffEqRODEAlgorithm
        @test RandomHeun() isa StochasticDiffEqRODEAlgorithm
        @test RandomTamedEM() isa StochasticDiffEqRODEAlgorithm
        @test BAOAB() isa StochasticDiffEqAlgorithm
        @test BAOAB(gamma = 2.0, scale_noise = false).gamma == 2.0
    end
end
