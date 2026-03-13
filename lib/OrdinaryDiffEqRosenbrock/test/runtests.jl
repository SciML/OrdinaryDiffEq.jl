using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "DAE Rosenbrock AD Tests" include("dae_rosenbrock_ad_tests.jl")
    @time @safetestset "Rosenbrock Convergence Tests" include("ode_rosenbrock_tests.jl")
    @time @safetestset "Jacobian Reuse Tests" include("jacobian_reuse_test.jl")
end

# Run QA tests (JET, Aqua, AllocCheck)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
