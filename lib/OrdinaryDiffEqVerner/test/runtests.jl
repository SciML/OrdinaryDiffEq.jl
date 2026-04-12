using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "RKV76IIa Tests" include("ode_verner_tests.jl")
end

# Run QA tests (JET, Aqua, AllocCheck)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
