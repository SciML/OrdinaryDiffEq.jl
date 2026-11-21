using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

if TEST_GROUP == "ALL" || TEST_GROUP == "Core"
    @time @safetestset "Module loads and constructors" begin
        using StochasticDiffEqImplicit
        using SciMLBase: SDEProblem, solve, successful_retcode
        using Test

        @test ImplicitEM() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitEulerHeun() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ImplicitRKMil() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test STrapezoid() isa ImplicitEM
        @test SImplicitMidpoint() isa ImplicitEM
        @test ISSEM() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test ISSEulerHeun() isa StochasticDiffEqNewtonAdaptiveAlgorithm
        @test SKenCarp() isa StochasticDiffEqNewtonAdaptiveAlgorithm

        f = (du, u, p, t) -> (du .= -u; nothing)
        g = (du, u, p, t) -> (du .= zero(eltype(u)); nothing)
        prob = SDEProblem(f, g, ones(2), (0.0, 1.0))
        sol = solve(prob, ISSEM(); dt = 0.1, adaptive = false)
        @test successful_retcode(sol)
    end
end

# Run QA tests (Aqua, JET) - skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "QA (Aqua and JET)" include("qa/qa.jl")
end
