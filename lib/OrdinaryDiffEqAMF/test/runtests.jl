using Test
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @testset "OrdinaryDiffEqAMF" begin
        include("test_pollu.jl")
        include("test_fd2d.jl")
        include("test_adjoint_fd2d.jl")
    end
end

# Run QA tests (AllocCheck) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if TEST_GROUP != "Core" && isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
