using Pkg
using SafeTestsets, Test
const LONGER_TESTS = false

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_odeinterface_env()
    Pkg.activate("odeinterface")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_ad_env()
    Pkg.activate("ad")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_modelingtoolkit_env()
    Pkg.activate("modelingtoolkit")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

#Start Test Script

@time begin
    # Detect sublibrary test groups.
    # GROUP can be a bare sublibrary name (Core test group) or
    # "{sublibrary}_{TEST_GROUP}" for any custom group (e.g., QA, GPU, etc.).
    # Sublibraries declare their groups in test/test_groups.toml.
    lib_dir = joinpath(dirname(@__DIR__), "lib")

    # Check if GROUP matches a sublibrary, possibly with a _SUFFIX for the test group.
    # Scan underscores right-to-left to find the longest matching sublibrary prefix.
    function _detect_sublibrary_group(group, lib_dir)
        isdir(joinpath(lib_dir, group)) && return (group, "Core")
        for i in length(group):-1:1
            if group[i] == '_' && isdir(joinpath(lib_dir, group[1:(i - 1)]))
                return (group[1:(i - 1)], group[(i + 1):end])
            end
        end
        return (group, "Core")
    end
    base_group, test_group = _detect_sublibrary_group(GROUP, lib_dir)

    if isdir(joinpath(lib_dir, base_group))
        Pkg.activate(joinpath(lib_dir, base_group))
        # On Julia < 1.11, the [sources] section in Project.toml is not supported.
        # Manually Pkg.develop local path dependencies so CI tests the PR branch code.
        if VERSION < v"1.11.0-DEV.0"
            toml = Pkg.TOML.parsefile(joinpath(lib_dir, base_group, "Project.toml"))
            if haskey(toml, "sources")
                for (dep_name, source_spec) in toml["sources"]
                    if source_spec isa Dict && haskey(source_spec, "path")
                        dep_path = normpath(joinpath(lib_dir, base_group, source_spec["path"]))
                        if isdir(dep_path)
                            @info "Developing local source dependency" dep_name dep_path
                            Pkg.develop(Pkg.PackageSpec(path = dep_path))
                        end
                    end
                end
            end
        end
        withenv("ODEDIFFEQ_TEST_GROUP" => test_group) do
            Pkg.test(base_group, julia_args = ["--check-bounds=auto", "--compiled-modules=yes", "--depwarn=yes"], force_latest_compatible_version = false, allow_reresolve = true)
        end
    elseif GROUP == "All" || GROUP == "InterfaceI" || GROUP == "Interface"
        @time @safetestset "Discrete Algorithm Tests" include("interface/discrete_algorithm_test.jl")
        # Skip on Julia LTS (oneunit(Type{Any}) not defined) and pre-release (stalls)
        # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/2979
        if VERSION >= v"1.11" && isempty(VERSION.prerelease)
            @time @safetestset "Null u0 Callbacks Tests" include("interface/null_u0_callbacks_test.jl")
        end
        @time @safetestset "Tstops Tests" include("interface/ode_tstops_tests.jl")
        @time @safetestset "Backwards Tests" include("interface/ode_backwards_test.jl")
        @time @safetestset "Initdt Tests" include("interface/ode_initdt_tests.jl")
        @time @safetestset "Linear Tests" include("interface/ode_twodimlinear_tests.jl")
        @time @safetestset "Inf Tests" include("interface/inf_handling.jl")
        @time @safetestset "saveat Tests" include("interface/ode_saveat_tests.jl")
        @time @safetestset "save_idxs Tests" include("interface/ode_saveidxs_tests.jl")
        @time @safetestset "Scalar Handling Tests" include("interface/scalar_handling_tests.jl")
        @time @safetestset "Static Array Tests" include("interface/static_array_tests.jl")
        @time @safetestset "u_modified Tests" include("interface/umodified_test.jl")
        @time @safetestset "Composite Algorithm Tests" include("interface/composite_algorithm_test.jl")
        @time @safetestset "Complex Tests" include("interface/complex_tests.jl")
        @time @safetestset "Ndim Complex Tests" include("interface/ode_ndim_complex_tests.jl")
        @time @safetestset "Number Type Tests" include("interface/ode_numbertype_tests.jl")
        @time @safetestset "Interpolation Output Type Tests" include("interface/interpolation_output_types.jl")
        @time @safetestset "Stiffness Detection Tests" include("interface/stiffness_detection_test.jl")
        @time @safetestset "Composite Interpolation Tests" include("interface/composite_interpolation.jl")
        @time @safetestset "Export tests" include("interface/export_tests.jl")
        @time @safetestset "Type Handling Tests" include("interface/type_handling.jl")
        @time @safetestset "Controller Tests" include("interface/controllers.jl")
        @time @safetestset "Inplace Interpolation Tests" include("interface/inplace_interpolation.jl")
        @time @safetestset "Algebraic Interpolation Tests" include("interface/algebraic_interpolation.jl")
        @time @safetestset "Interpolation and Cache Stripping Tests" include("interface/ode_strip_test.jl")
        @time @safetestset "Aliasing Tests" include("interface/aliasing_tests.jl")
        @time @safetestset "Solution Memory Release" include("interface/solution_memory_tests.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceII" || GROUP == "Interface")
        #@time @safetestset "No Recompile Tests" include("interface/norecompile.jl") # doesn't work on CI?
        @time @safetestset "AutoSparse Detection Tests" include("interface/autosparse_detection_tests.jl")
        @time @safetestset "Enum Tests" include("interface/enums.jl")
        @time @safetestset "Get du Tests" include("interface/get_du.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceIII" || GROUP == "Interface")
        @time @safetestset "Derivative Utilities Tests" include("interface/utility_tests.jl")
        @time @safetestset "stats Tests" include("interface/stats_tests.jl")
        @time @safetestset "No Index Tests" include("interface/noindex_tests.jl")
        @time @safetestset "Events + DAE addsteps Tests" include("interface/event_dae_addsteps.jl")
        @time @safetestset "Units Tests" include("interface/units_tests.jl")
        @time @safetestset "DEVerbosity Tests" include("interface/verbosity.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceIV" || GROUP == "Interface")
        @time @safetestset "Ambiguity Tests" include("interface/ambiguity_tests.jl")
        @time @safetestset "Precision Mixing Tests" include("interface/precision_mixing.jl")
        @time @safetestset "Sized Matrix Tests" include("interface/sized_matrix_tests.jl")
        @time @safetestset "Second Order with First Order Solver Tests" include("interface/second_order_with_first_order_solvers.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceV" || GROUP == "Interface")
        @time @safetestset "Interpolation Derivative Error Tests" include("interface/interpolation_derivative_error_tests.jl")
        @time @safetestset "GPU AutoDiff Interface Tests" include("interface/gpu_autodiff_interface_tests.jl")
    end

    if !is_APPVEYOR &&
            (GROUP == "All" || GROUP == "Integrators_I" || GROUP == "Integrators")
        @time @safetestset "Reinit Tests" include("integrators/reinit_test.jl")
        @time @safetestset "Events Tests" include("integrators/ode_event_tests.jl")
        @time @safetestset "Alg Events Tests" include("integrators/alg_events_tests.jl")
        @time @safetestset "Discrete Callback Dual Tests" include("integrators/discrete_callback_dual_test.jl")
        @time @safetestset "Callback Allocation Tests" include("integrators/callback_allocation_tests.jl")
        @time @safetestset "Iterator Tests" include("integrators/iterator_tests.jl")
        @time @safetestset "Integrator Interface Tests" include("integrators/integrator_interface_tests.jl")
        @time @safetestset "Error Check Tests" include("integrators/check_error.jl")
        @time @safetestset "Event Detection Tests" include("integrators/event_detection_tests.jl")
        @time @safetestset "Event Repetition Detection Tests" include("integrators/event_repeat_tests.jl")
        @time @safetestset "Step Limiter Tests" include("integrators/step_limiter_test.jl")
    end

    if !is_APPVEYOR &&
            (GROUP == "All" || GROUP == "Integrators_II" || GROUP == "Integrators")
        @time @safetestset "Integrator RNG Tests" include("integrators/integrator_rng_tests.jl")
        @time @safetestset "Reverse Directioned Event Tests" include("integrators/rev_events_tests.jl")
        @time @safetestset "Differentiation Direction Tests" include("integrators/diffdir_tests.jl")
        @time @safetestset "Resize Tests" include("integrators/resize_tests.jl")
        @time @safetestset "Cache Tests" include("integrators/ode_cache_tests.jl")
        @time @safetestset "Add Steps Tests" include("integrators/ode_add_steps_tests.jl")
        @time @safetestset "IMEX Split Function Tests" include("integrators/split_ode_tests.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "Regression_I" || GROUP == "Regression")
        @time @safetestset "Dense Tests" include("regression/ode_dense_tests.jl")
        @time @safetestset "Special Interp Tests" include("regression/special_interps.jl")
        @time @safetestset "Inplace Tests" include("regression/ode_inplace_tests.jl")
        @time @safetestset "Adaptive Tests" include("regression/ode_adaptive_tests.jl")
        @time @safetestset "Hard DAE Tests" include("regression/hard_dae.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "Regression_II" || GROUP == "Regression")
        @time @safetestset "PSOS Energy Conservation Tests" include("regression/psos_and_energy_conservation.jl")
        @time @safetestset "Unrolled Tests" include("regression/ode_unrolled_comparison_tests.jl")
        @time @safetestset "IIP vs OOP Tests" include("regression/iipvsoop_tests.jl")
        @time @safetestset "Inference Tests" include("regression/inference.jl")
    end

    if !is_APPVEYOR && GROUP == "AlgConvergence_I"
        @time @safetestset "Non-autonomous Convergence Tests" include("algconvergence/non-autonomous_convergence_tests.jl")
    end

    if !is_APPVEYOR && GROUP == "AlgConvergence_III"
        @time @safetestset "Split Methods Tests" include("algconvergence/split_methods_tests.jl")
    end

    # ModelingToolkit tests moved to OrdinaryDiffEqDifferentiation and
    # OrdinaryDiffEqNonlinearSolve subpackage test groups (ModelingToolkit group).

    if !is_APPVEYOR && GROUP == "Downstream"
        activate_downstream_env()
        @time @safetestset "DelayDiffEq Tests" include("downstream/delaydiffeq.jl")
        @time @safetestset "Measurements Tests" include("downstream/measurements.jl")
        @time @safetestset "Time derivative Tests" include("downstream/time_derivative_test.jl")
    end

    # AD tests - Enzyme/Zygote only on Julia <= 1.11 (see https://github.com/EnzymeAD/Enzyme.jl/issues/2699)
    # Mooncake works on all Julia versions
    if !is_APPVEYOR && GROUP == "AD"
        activate_ad_env()
        @time @safetestset "AD Tests" include("ad/ad_tests.jl")
        @time @safetestset "Autodiff Events Tests" include("ad/autodiff_events.jl")
        @time @safetestset "Discrete Adjoint Tests" include("ad/discrete_adjoints.jl")
    end

    # Don't run ODEInterface tests on prerelease
    if !is_APPVEYOR && GROUP == "ODEInterfaceRegression" && isempty(VERSION.prerelease)
        activate_odeinterface_env()
        @time @safetestset "Init dt vs dorpri tests" include("odeinterface/init_dt_vs_dopri_tests.jl")
        @time @safetestset "ODEInterface Regression Tests" include("odeinterface/odeinterface_regression.jl")
    end

    if !is_APPVEYOR && GROUP == "QA"
        @time @safetestset "Quality Assurance Tests" include("qa/qa_tests.jl")
    end
end # @time
