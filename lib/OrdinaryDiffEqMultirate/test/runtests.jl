using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "MREEF Tests" include("mreef_tests.jl")
end

# Run QA tests (AllocCheck) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
