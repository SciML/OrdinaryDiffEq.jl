using SciMLTesting
using SafeTestsets

const TEST_GROUP = get(ENV, "GROUP", "ALL")

function activate_qa_env()
    return activate_group_env(
        joinpath(@__DIR__, "qa");
        parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")]
    )
end

if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "SDC Tableau Tests" include("sdc_tableau_tests.jl")
    @time @safetestset "SDC Convergence Tests" include("sdc_convergence_tests.jl")
    @time @safetestset "SDC Stiff Tests" include("sdc_stiff_tests.jl")
    @time @safetestset "SDC Diagonal Sweeper Tests" include("sdc_diagonal_sweeper_tests.jl")
    @time @safetestset "SDC Adaptive Tests" include("sdc_adaptive_tests.jl")
end

# Allocation tests must run before JET: JET's static analysis invalidates
# compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "Quality Assurance" include("qa/qa.jl")
    @time @safetestset "JET" include("qa/jet.jl")
end
