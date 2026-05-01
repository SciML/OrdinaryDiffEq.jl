using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

function activate_threaded_env()
    Pkg.activate(joinpath(@__DIR__, "threaded"))
    return Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "ABM Convergence Tests" include("abm_convergence_tests.jl")
    @time @safetestset "Adams Variable Coefficients Tests" include("adams_tests.jl")
end

# Threaded tests require Polyester.jl (for FastBroadcast.Threaded() support).
# Runs only on explicit `ODEDIFFEQ_TEST_GROUP=Threaded` — the SublibraryCI
# matrix (lib/OrdinaryDiffEqAdamsBashforthMoulton/test/test_groups.toml)
# schedules this as its own Julia 1 / 2-thread job. Matches the GPU-group
# convention in sibling sublibs (LowStorageRK, Rosenbrock, BDF).
if TEST_GROUP == "Threaded"
    activate_threaded_env()
    @time @safetestset "ABM Threaded Convergence Tests" include("threaded/abm_threaded_convergence_tests.jl")
    @time @safetestset "Regression test for threading versions vs non threading versions" include("threaded/regression_test_threading.jl")
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
