using DiffEqDevTools
using Pkg
using SafeTestsets
using Test

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# Run QA tests (Aqua) — skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Aqua" include("qa/qa.jl")
end

if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    # write your own tests here
    @time @safetestset "Benchmark Tests" include("benchmark_tests.jl")
    @time @safetestset "ODE AppxTrue Tests" include("ode_appxtrue_tests.jl")
    @time @safetestset "Analyticless Convergence Tests" include("analyticless_convergence_tests.jl")
    @time @safetestset "Analyticless Stochastic WP" include("analyticless_stochastic_wp.jl")
    @time @safetestset "Stability Region Tests" include("stability_region_test.jl")
    @time @safetestset "Plot Recipes" include("plotrecipes_tests.jl")
    @time @safetestset "Plot Recipes (Nonlinearsolve WP-diagrams)" include("nonlinearsolve_wpdiagram_tests.jl")
end
