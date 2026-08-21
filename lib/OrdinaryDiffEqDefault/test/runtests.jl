using Pkg
using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_gpu_env()
    Pkg.activate(joinpath(@__DIR__, "gpu"))
    return Pkg.instantiate()
end

# Run GPU tests
if TEST_GROUP == "GPU"
    activate_gpu_env()
    @time @safetestset "Autoswitch GPU" include("gpu/autoswitch.jl")
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "SciMLBase reexport" begin
        # docs/src/api/reexports.md defines this surface; test/qa/qa_tests.jl checks the
        # full list repo-wide. Here, spot-check that the common interface is usable and
        # that solver-author API stayed behind the `SciMLBase.` qualifier.
        using OrdinaryDiffEqDefault, Test
        exported = (
            :ODEProblem, :ODEFunction, :SplitODEProblem, :solve, :init, :step!,
            :remake, :ReturnCode, :CallbackSet, :ContinuousCallback, :terminate!,
            :u_modified!, :add_tstop!, :get_du, :EnsembleProblem,
        )
        @test all(Base.isexported.(Ref(OrdinaryDiffEqDefault), exported))
        internal = (
            :build_solution, :isinplace, :has_jac, :AbstractODEProblem,
            :StandardODEProblem, :UJacobianWrapper,
        )
        @test !any(Base.isexported.(Ref(OrdinaryDiffEqDefault), internal))
    end
    @time @safetestset "Default Solver Tests" include("default_solver_tests.jl")
end

# Run QA tests (JET, Aqua)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
