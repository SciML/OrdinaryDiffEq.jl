using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Low Order ERK Convergence Tests" include("low_order_erk_convergence_tests.jl")
    @time @safetestset "OwrenZen Tests" include("owrenzen_tests.jl")
    @time @safetestset "Euler SSP Tests" include("euler_ssp.jl")
end

# Run QA tests (AllocCheck, JET, Aqua)
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
