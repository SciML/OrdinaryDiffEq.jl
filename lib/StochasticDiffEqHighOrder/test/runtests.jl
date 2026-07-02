using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

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

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
