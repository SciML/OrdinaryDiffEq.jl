using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")
const ORIGINAL_PROJECT = dirname(Base.active_project())

function activate_gpu_env()
    Pkg.activate(joinpath(@__DIR__, "gpu"))
    return Pkg.instantiate()
end

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.instantiate()
end

# Run GPU tests
if TEST_GROUP == "GPU"
    activate_gpu_env()
    @time @safetestset "BDF Solvers GPU" include("gpu/bdf_solvers.jl")
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "DAE Convergence Tests" include("dae_convergence_tests.jl")
    @time @safetestset "DAE AD Tests" include("dae_ad_tests.jl")
    @time @safetestset "DAE Event Tests" include("dae_event.jl")
    @time @safetestset "DAE derivative_discontinuity! Tests" include("dae_derivative_discontinuity_tests.jl")
    @time @safetestset "DAE Initialization Tests" include("dae_initialization_tests.jl")

    @time @safetestset "BDF Inference Tests" include("inference_tests.jl")
    @time @safetestset "BDF Convergence Tests" include("bdf_convergence_tests.jl")
    @time @safetestset "BDF Regression Tests" include("bdf_regression_tests.jl")
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
    Pkg.activate(ORIGINAL_PROJECT)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
