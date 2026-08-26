using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

@time @safetestset "SciMLBase reexport" begin
    using OrdinaryDiffEqPDIRK, Test
    exported = (
        :ODEProblem, :ODEFunction, :SplitODEProblem, :solve, :init, :step!,
        :remake, :ReturnCode, :CallbackSet, :ContinuousCallback, :terminate!,
        :u_modified!, :add_tstop!, :get_du, :EnsembleProblem,
    )
    @test all(Base.isexported.(Ref(OrdinaryDiffEqPDIRK), exported))
    internal = (
        :build_solution, :isinplace, :has_jac, :AbstractODEProblem,
        :StandardODEProblem, :UJacobianWrapper, :LinearProblem,
        :ConvexOptimizationProblem,
    )
    @test !any(Base.isexported.(Ref(OrdinaryDiffEqPDIRK), internal))
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
