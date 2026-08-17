using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

@time @safetestset "Discontinuity restart" include("imex_discontinuity_restart_tests.jl")

@time @safetestset "SciMLBase reexport" begin
    using OrdinaryDiffEqIMEXMultistep, Test
    expected = (
        :ODEProblem, :ODEFunction, :solve, :init, :solve!, :step!, :remake, :reinit!,
        :ReturnCode, :ContinuousCallback, :DiscreteCallback, :VectorContinuousCallback,
        :CallbackSet, :terminate!, :add_tstop!, :derivative_discontinuity!,
        :set_proposed_dt!, :successful_retcode, :ODEAliasSpecifier,
        :SplitODEProblem, :SplitFunction,
    )
    @test all(Base.isexported.(Ref(OrdinaryDiffEqIMEXMultistep), expected))
    @test !Base.isexported(OrdinaryDiffEqIMEXMultistep, :EnsembleProblem)
    @test !Base.isexported(OrdinaryDiffEqIMEXMultistep, :get_du)
    @test !Base.isexported(OrdinaryDiffEqIMEXMultistep, :u_modified!)
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
