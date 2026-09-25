using Pkg
using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_gpu_env()
    Pkg.activate(joinpath(@__DIR__, "gpu"))
    return Pkg.instantiate()
end

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run GPU tests
if TEST_GROUP == "GPU"
    activate_gpu_env()
    @time @safetestset "Simple GPU" include("gpu/simple_gpu.jl")
    @time @safetestset "Hermite Interpolation GPU" include("gpu/hermite_test.jl")
end

# Functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Generic algorithm trait contract" include("algorithm_interface_tests.jl")
    @time @safetestset "change_t_via_interpolation!" include("change_t_via_interpolation_tests.jl")
    @time @safetestset "SciMLBase constructor bindings" begin
        using OrdinaryDiffEqCore: DynamicalODEProblem, ODEFunction
        using SciMLBase, Test

        @test ODEFunction === SciMLBase.ODEFunction
        @test DynamicalODEProblem === SciMLBase.DynamicalODEProblem
    end
    @time @safetestset "Developer Time Queue API" include("developer_time_queue_api_tests.jl")
    @time @safetestset "Developer Codegen API" include("developer_codegen_api_tests.jl")
    @time @safetestset "Sparse isdiag Performance" include("sparse_isdiag_tests.jl")
    @time @safetestset "Algebraic Vars Detection" include("algebraic_vars_detection_tests.jl")
    @time @safetestset "Interpolation Search Hint" include("interpolation_hint_tests.jl")
    @time @safetestset "Bool Equal Coercion" include("bool_equal_tests.jl")
    @time @safetestset "dtmin Direction" include("dtmin_direction_tests.jl")
    @time @safetestset "Instability Diagnostics" include("instability_diagnostics_tests.jl")
    @time @safetestset "Enzyme Interpolation" include("enzyme_interpolation_tests.jl")
end

# Run QA tests LAST. `JET.test_package` re-evaluates this package's source into a
# virtual module, so every method the package defines on a generic function owned by
# another module is replaced by a copy bound to a module with no package extensions
# loaded. Anything that runs afterwards in the same process then exercises those
# copies instead of the real methods. `activate_qa_env()` also leaves the QA
# environment active, so the groups above must resolve before it runs.
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Core Infrastructure AllocCheck Tests" include("qa/alloccheck.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
