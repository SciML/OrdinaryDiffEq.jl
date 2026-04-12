using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_threaded_env()
    Pkg.activate(joinpath(@__DIR__, "threaded"))
    return Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "ABM Convergence Tests" include("abm_convergence_tests.jl")
    @time @safetestset "Adams Variable Coefficients Tests" include("adams_tests.jl")
end

# Threaded tests require Polyester.jl (for FastBroadcast.Threaded() support)
if TEST_GROUP == "Threaded" || TEST_GROUP == "ALL"
    activate_threaded_env()
    @time @safetestset "ABM Threaded Convergence Tests" include("threaded/abm_threaded_convergence_tests.jl")
    @time @safetestset "Regression test for threading versions vs non threading versions" include("threaded/regression_test_threading.jl")
end

# Run QA tests (JET, Aqua)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
