using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

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

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
