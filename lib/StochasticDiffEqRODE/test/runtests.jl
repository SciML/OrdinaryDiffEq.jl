using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

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

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
