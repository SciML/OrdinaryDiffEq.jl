using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")
const ORIGINAL_PROJECT = dirname(Base.active_project())

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Extrapolation Tests" include("ode_extrapolation_tests.jl")
end

if TEST_GROUP == "Multithreading"
    @time @safetestset "Multithreaded Extrapolation Tests" include("multithreading/ode_extrapolation_tests.jl")
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
    Pkg.activate(ORIGINAL_PROJECT)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
