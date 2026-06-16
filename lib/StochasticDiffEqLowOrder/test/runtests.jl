using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads" begin
        using StochasticDiffEqLowOrder
        using Test

        @test isdefined(StochasticDiffEqLowOrder, :EM)
        @test isdefined(StochasticDiffEqLowOrder, :EulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEM)
        @test isdefined(StochasticDiffEqLowOrder, :LambaEulerHeun)
        @test isdefined(StochasticDiffEqLowOrder, :SimplifiedEM)
        @test isdefined(StochasticDiffEqLowOrder, :SplitEM)
        @test isdefined(StochasticDiffEqLowOrder, :RKMil)
        @test isdefined(StochasticDiffEqLowOrder, :RKMilCommute)
        @test isdefined(StochasticDiffEqLowOrder, :PCEuler)
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
