using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end

@time @safetestset "Convergence" include("sdirk_convergence_tests.jl")
@time @safetestset "DAE tests" include("dae_esdirk_test.jl")
