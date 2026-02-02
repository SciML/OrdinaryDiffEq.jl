using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Convergence Tests" begin
        include("convergence_tests.jl")
    end
    @time @safetestset "Adaptivity Tests" begin
        include("adaptivity_tests.jl")
    end
end

# Run QA tests (JET)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
