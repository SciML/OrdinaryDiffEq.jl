using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

@time @safetestset "SciMLBase reexport" begin
    using OrdinaryDiffEqPDIRK, Test
    expected = (
        :ODEProblem, :ODEFunction, :solve, :init, :solve!, :step!, :remake, :reinit!,
        :ReturnCode, :ContinuousCallback, :DiscreteCallback, :VectorContinuousCallback,
        :CallbackSet, :terminate!, :add_tstop!, :derivative_discontinuity!,
        :set_proposed_dt!, :successful_retcode, :ODEAliasSpecifier,
    )
    @test all(Base.isexported.(Ref(OrdinaryDiffEqPDIRK), expected))
    @test !Base.isexported(OrdinaryDiffEqPDIRK, :EnsembleProblem)
    @test !Base.isexported(OrdinaryDiffEqPDIRK, :get_du)
    @test !Base.isexported(OrdinaryDiffEqPDIRK, :u_modified!)
end

# Run QA tests (AllocCheck, JET) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end

@time @safetestset "Convergence Tests" include("pdirk_convergence_tests.jl")
@time @safetestset "nlsolve! Arguments" include("nlsolve_argument_tests.jl")
