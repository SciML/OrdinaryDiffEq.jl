using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP != "QA"
    @time @safetestset "MREEF Tests" include("mreef_tests.jl")
end
