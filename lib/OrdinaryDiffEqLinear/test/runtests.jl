using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Linear Methods Tests" include("linear_method_tests.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL"
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
