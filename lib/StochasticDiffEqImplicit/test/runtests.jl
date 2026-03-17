using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqImplicit
        using Test

        @test ImplicitEM() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitEulerHeun() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitRKMil() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test STrapezoid() isa ImplicitEM
        @test SImplicitMidpoint() isa ImplicitEM
        @test ISSEM() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ISSEulerHeun() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test SKenCarp() isa StochasticDiffEqNewtonAdaptiveAlgorithm
    end
end
