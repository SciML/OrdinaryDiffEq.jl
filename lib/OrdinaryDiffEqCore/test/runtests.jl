using Pkg
using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

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

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Core Infrastructure AllocCheck Tests" include("qa/alloccheck.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end

# Functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Sparse isdiag Performance" include("sparse_isdiag_tests.jl")
    @time @safetestset "Algebraic Vars Detection" include("algebraic_vars_detection_tests.jl")
    @time @safetestset "Discontinuity Detection" include("disco_tests.jl")
    @time @safetestset "Interpolation Search Hint" include("interpolation_hint_tests.jl")
end
