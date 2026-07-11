using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

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

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
