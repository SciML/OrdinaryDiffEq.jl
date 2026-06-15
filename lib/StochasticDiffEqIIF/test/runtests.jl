using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqIIF
        using Test

        @test IIF1M() isa IIF1M
        @test IIF2M() isa IIF2M
        @test IIF1Mil() isa IIF1Mil
        @test StochasticDiffEqIIF.alg_order(IIF1M()) == 1 // 2
        @test StochasticDiffEqIIF.alg_order(IIF2M()) == 1 // 2
        @test StochasticDiffEqIIF.alg_order(IIF1Mil()) == 1 // 1
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
