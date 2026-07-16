using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run QA tests (AllocCheck, JET, Aqua) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
end

@time @safetestset "Tableau consistency" include("tableau_consistency_tests.jl")
@time @safetestset "Stage predictors" include("predictor_tests.jl")
@time @safetestset "Convergence" include("sdirk_convergence_tests.jl")
@time @safetestset "DAE tests" include("dae_esdirk_test.jl")
