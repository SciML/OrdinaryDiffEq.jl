using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run QA tests (JET, Aqua, AllocCheck)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
