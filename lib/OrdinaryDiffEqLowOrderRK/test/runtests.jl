using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Low Order ERK Convergence Tests" include("low_order_erk_convergence_tests.jl")
    @time @safetestset "OwrenZen Tests" include("owrenzen_tests.jl")
    @time @safetestset "Euler SSP Tests" include("euler_ssp.jl")
end

# Run QA tests (JET, Aqua, AllocCheck)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
