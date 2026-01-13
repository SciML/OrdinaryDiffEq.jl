using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL"
    @time @safetestset "JET Tests" include("jet.jl")
end
