using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Convergence Tests" begin
        include("convergence_tests.jl")
    end
    @time @safetestset "Adaptivity Tests" begin
        include("adaptivity_tests.jl")
    end
end

# Run QA tests (AllocCheck, JET) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
    @time @safetestset "JET Tests" include("jet.jl")
end
