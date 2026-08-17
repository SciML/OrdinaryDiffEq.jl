using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "SciMLBase reexport" begin
        using OrdinaryDiffEqFunctionMap, Test
        expected = (
            :ODEProblem, :ODEFunction, :solve, :init, :solve!, :step!, :remake, :reinit!,
            :ReturnCode, :ContinuousCallback, :DiscreteCallback, :VectorContinuousCallback,
            :CallbackSet, :terminate!, :add_tstop!, :derivative_discontinuity!,
            :set_proposed_dt!, :successful_retcode, :ODEAliasSpecifier,
            :DiscreteProblem, :DiscreteFunction,
        )
        @test all(Base.isexported.(Ref(OrdinaryDiffEqFunctionMap), expected))
        @test !Base.isexported(OrdinaryDiffEqFunctionMap, :EnsembleProblem)
        @test !Base.isexported(OrdinaryDiffEqFunctionMap, :get_du)
        @test !Base.isexported(OrdinaryDiffEqFunctionMap, :u_modified!)
    end
    @time @safetestset "DiscreteProblem Defaults" include("discrete_problem_defaults.jl")
    @time @safetestset "Discrete Algorithm Tests" include("discrete_algorithm_test.jl")
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
