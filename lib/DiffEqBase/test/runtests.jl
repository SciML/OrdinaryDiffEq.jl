using Pkg
using SafeTestsets
using Test

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "Core")

function activate_downstream_env()
    Pkg.activate(joinpath(@__DIR__, "downstream"))
    return Pkg.instantiate()
end

function activate_static_env()
    Pkg.activate(joinpath(@__DIR__, "static"))
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_modelingtoolkit_env()
    Pkg.activate(joinpath(@__DIR__, "modelingtoolkit"))
    return Pkg.instantiate()
end

function activate_sundials_env()
    Pkg.activate(joinpath(@__DIR__, "sundials"))
    return Pkg.instantiate()
end

@time begin
    # Core tests — basic DiffEqBase functionality, no downstream solvers needed
    if TEST_GROUP ∉ (
            "QA", "Static", "Downstream", "Downstream2",
            "ModelingToolkit", "Sundials",
        )
        @time @safetestset "Callbacks" include("callbacks.jl")
        @time @safetestset "Plot Vars" include("plot_vars.jl")
        @time @safetestset "Problem Creation Tests" include("problem_creation_tests.jl")
        @time @safetestset "Export tests" include("export_tests.jl")
        @time @safetestset "Remake tests" include("remake_tests.jl")
        @time @safetestset "High Level solve Interface" include("high_level_solve.jl")
        @time @safetestset "DiffEqFunction tests" include("diffeqfunction_tests.jl")
        @time @safetestset "Internal Euler" include("internal_euler_test.jl")
        @time @safetestset "Norm" include("norm.jl")
        @time @safetestset "Utils" include("utils.jl")
        @time @safetestset "ForwardDiff Dual Detection" include("forwarddiff_dual_detection.jl")
        @time @safetestset "ODE default norm" include("ode_default_norm.jl")
        @time @safetestset "DynamicQuantities extension" include("dynamicquantities_ext.jl")
        @time @safetestset "ODE default unstable check" include("ode_default_unstable_check.jl")
        @time @safetestset "Problem Kwargs Merging" include("problem_kwargs_merging.jl")
        @time @safetestset "Verbose Inference" include("verbose_inference.jl")
    end

    # QA tests — Aqua quality checks
    if TEST_GROUP ∉ (
            "Core", "Static", "Downstream", "Downstream2",
            "ModelingToolkit", "Sundials",
        ) && isempty(VERSION.prerelease)
        @time @safetestset "Aqua" include("aqua.jl")
    end

    # Static analysis tests — allocation checks with ComponentArrays
    if TEST_GROUP == "Static" || TEST_GROUP == "ALL"
        activate_static_env()
        @time @safetestset "Static Checks" include("static/static_checks.jl")
    end

    # Downstream tests — OrdinaryDiffEq integration (lightweight deps)
    if TEST_GROUP == "Downstream" || TEST_GROUP == "ALL"
        activate_downstream_env()
        @time @safetestset "Kwarg Warnings" include("downstream/kwarg_warn.jl")
        @time @safetestset "Solve Error Handling" include("downstream/solve_error_handling.jl")
        @time @safetestset "Dual Detection Solution" include("downstream/dual_detection_solution.jl")
        @time @safetestset "Null Parameters" include("downstream/null_params_test.jl")
        @time @safetestset "Ensemble Simulations" include("downstream/ensemble.jl")
        @time @safetestset "Ensemble Analysis" include("downstream/ensemble_analysis.jl")
        @time @safetestset "Ensemble Thread Safety" include("downstream/ensemble_thread_safety.jl")
        @time @safetestset "Inference Tests" include("downstream/inference.jl")
        @time @safetestset "Table Inference Tests" include("downstream/tables.jl")
        @time @safetestset "Default linsolve with structure" include("downstream/default_linsolve_structure.jl")
        @time @safetestset "Callback Merging Tests" include("downstream/callback_merging.jl")
        @time @safetestset "Callback Detection Tests" include("downstream/callback_detection.jl")
        @time @safetestset "LabelledArrays Tests" include("downstream/labelledarrays.jl")
        @time @safetestset "GTPSA Tests" include("downstream/gtpsa.jl")
        @time @safetestset "SubArray Support" include("downstream/subarray_support.jl")
        @time @safetestset "Unitful" include("downstream/unitful.jl")
        @time @safetestset "FlexUnits" include("downstream/flexunits.jl")
    end

    # Downstream2 tests — additional OrdinaryDiffEq integration tests
    if TEST_GROUP == "Downstream2" || TEST_GROUP == "ALL"
        if TEST_GROUP != "Downstream"
            activate_downstream_env()
        end
        @time @safetestset "Prob Kwargs" include("downstream/prob_kwargs.jl")
        @time @safetestset "Unwrapping" include("downstream/unwrapping.jl")
        @time @safetestset "Callback BigFloats" include("downstream/bigfloat_events.jl")
        @time @safetestset "DE stats" include("downstream/stats_tests.jl")
        @time @safetestset "Community Callback Tests" include("downstream/community_callback_tests.jl")
        @time @testset "Distributed Ensemble Tests" include("downstream/distributed_ensemble.jl")
    end

    # ModelingToolkit tests — heavy MTK dependency, separate environment
    if TEST_GROUP == "ModelingToolkit" && isempty(VERSION.prerelease)
        activate_modelingtoolkit_env()
        @time @safetestset "Null DE Handling" include("modelingtoolkit/null_de.jl")
    end

    # Sundials tests — only run when DiffEqBase itself changes
    if TEST_GROUP == "Sundials"
        activate_sundials_env()
        @time @safetestset "Sundials Error Handling" include("sundials/sundials_error_handling.jl")
    end
end
