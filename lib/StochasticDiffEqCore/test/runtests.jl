using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and types" begin
        using StochasticDiffEqCore
        using Test

        @testset "Module loads" begin
            @test isdefined(StochasticDiffEqCore, :SDEIntegrator)
            @test isdefined(StochasticDiffEqCore, :StochasticDiffEqAlgorithm)
            @test isdefined(StochasticDiffEqCore, :StochasticDiffEqAdaptiveAlgorithm)
            @test isdefined(StochasticDiffEqCore, :StochasticCompositeAlgorithm)
            @test isdefined(StochasticDiffEqCore, :IICommutative)
            @test isdefined(StochasticDiffEqCore, :IILevyArea)
            @test isdefined(StochasticDiffEqCore, :AutoSwitch)
        end

        @testset "Abstract type hierarchy" begin
            @test StochasticDiffEqNewtonAdaptiveAlgorithm <: StochasticDiffEqAdaptiveAlgorithm
            @test StochasticDiffEqNewtonAlgorithm <: StochasticDiffEqAlgorithm
            @test StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm
            @test StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm
        end

        @testset "Iterated integrals" begin
            @test IICommutative() isa IteratedIntegralApprox
            @test IILevyArea() isa IteratedIntegralApprox
        end
    end
end
