using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run QA tests (JET, Aqua)
if TEST_GROUP != "Core" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
