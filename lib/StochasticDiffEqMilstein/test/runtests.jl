using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

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

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
