using Pkg
using Test
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @testset "OrdinaryDiffEqAMF" begin
        include("test_pollu.jl")
        include("test_fd2d.jl")
        include("test_adjoint_fd2d.jl")
    end
end

# Run QA tests (AllocCheck) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
end
