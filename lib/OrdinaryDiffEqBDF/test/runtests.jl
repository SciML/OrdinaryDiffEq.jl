using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "DAE Convergence Tests" include("dae_convergence_tests.jl")
    @time @safetestset "DAE AD Tests" include("dae_ad_tests.jl")
    @time @safetestset "DAE Event Tests" include("dae_event.jl")
    @time @safetestset "DAE u_modified! Tests" include("dae_u_modified_tests.jl")
    @time @safetestset "DAE Initialization Tests" include("dae_initialization_tests.jl")

    @time @safetestset "BDF Inference Tests" include("inference_tests.jl")
    @time @safetestset "BDF Convergence Tests" include("bdf_convergence_tests.jl")
    @time @safetestset "BDF Regression Tests" include("bdf_regression_tests.jl")
    @time @safetestset "StiffInitDt Tests" include("stiff_initdt_tests.jl")
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
