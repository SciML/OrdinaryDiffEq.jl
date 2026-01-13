using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "ABM Convergence Tests" include("abm_convergence_tests.jl")
    @time @safetestset "Adams Variable Coefficients Tests" include("adams_tests.jl")
    @time @safetestset "Regression test for threading versions vs non threading versions" include("regression_test_threading.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL"
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
