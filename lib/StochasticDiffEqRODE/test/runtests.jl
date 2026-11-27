using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
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

        for alg in (RandomEM(), RandomHeun(), RandomTamedEM())
            @test StochasticDiffEqRODE.alg_order(alg) == 1 // 1
        end
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
