using Pkg
using SafeTestsets
using Test

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = (Sys.iswindows() && haskey(ENV, "APPVEYOR"))

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_static_env()
    Pkg.activate("static")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_gpu_env()
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

@time begin
    if GROUP == "All" || GROUP == "Core"
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
    end

    if !is_APPVEYOR && GROUP == "Downstream"
        activate_downstream_env()
        @time @safetestset "Kwarg Warnings" include("downstream/kwarg_warn.jl")
        @time @safetestset "Solve Error Handling" include("downstream/solve_error_handling.jl")
        @time @safetestset "Null DE Handling" include("downstream/null_de.jl")
        @time @safetestset "StaticArrays + AD" include("downstream/static_arrays_ad.jl")
        @time @safetestset "Unitful" include("downstream/unitful.jl")
        @time @safetestset "FlexUnits" include("downstream/flexunits.jl")
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
    end

    if !is_APPVEYOR && GROUP == "Static"
        activate_static_env()
        @time @safetestset "Static Checks" include("static/static_checks.jl")
    end

    if !is_APPVEYOR && GROUP == "Downstream2"
        activate_downstream_env()
        @time @safetestset "Prob Kwargs" include("downstream/prob_kwargs.jl")
        @time @safetestset "Unwrapping" include("downstream/unwrapping.jl")
        @time @safetestset "Callback BigFloats" include("downstream/bigfloat_events.jl")
        @time @safetestset "DE stats" include("downstream/stats_tests.jl")
        isempty(VERSION.prerelease) &&
            @time @safetestset "Ensemble AD Tests" include("downstream/ensemble_ad.jl")
        @time @safetestset "Community Callback Tests" include("downstream/community_callback_tests.jl")
        @time @safetestset "AD via ode with complex numbers" include("downstream/complex_number_ad.jl")
        @time @safetestset "Callback AD Tests" include("downstream/callback_ad.jl")
        @time @testset "Distributed Ensemble Tests" include("downstream/distributed_ensemble.jl")
    end

    if !is_APPVEYOR && GROUP == "GPU"
        activate_gpu_env()
        @time @safetestset "Simple GPU" include("gpu/simple_gpu.jl")
    end
end
