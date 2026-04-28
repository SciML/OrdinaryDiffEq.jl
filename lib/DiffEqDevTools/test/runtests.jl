using DiffEqDevTools
using Pkg
using Test

# SublibraryCI sets GROUP to "<pkg>_<group>" (e.g. "DiffEqDevTools_QA").
# CI.yml-style workflows set ODEDIFFEQ_TEST_GROUP. Default to ALL when neither
# is set so local `Pkg.test()` runs the full suite.
const RAW_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", get(ENV, "GROUP", "ALL"))
const TEST_GROUP = endswith(RAW_GROUP, "_QA") ? "QA" :
    (RAW_GROUP == "DiffEqDevTools" ? "Core" : RAW_GROUP)

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# Run QA tests (Aqua) — skip on pre-release Julia
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @testset "Aqua" begin
        include("qa/qa.jl")
    end
end

if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    # write your own tests here
    @time @testset "Benchmark Tests" begin
        include("benchmark_tests.jl")
    end
    @time @testset "ODE AppxTrue Tests" begin
        include("ode_appxtrue_tests.jl")
    end
    @time @testset "Analyticless Convergence Tests" begin
        include("analyticless_convergence_tests.jl")
    end
    @time @testset "Analyticless Stochastic WP" begin
        include("analyticless_stochastic_wp.jl")
    end
    @time @testset "Stability Region Tests" begin
        include("stability_region_test.jl")
    end
    @time @testset "Plot Recipes" begin
        include("plotrecipes_tests.jl")
    end
    @time @testset "Plot Recipes (Nonlinearsolve WP-diagrams)" begin
        include("nonlinearsolve_wpdiagram_tests.jl")
    end
end
