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

function activate_gpu_env()
    Pkg.activate("gpu")
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
    # Handle sublibrary QA groups (e.g., OrdinaryDiffEqBDF_QA)
    is_qa_group = endswith(GROUP, "_QA")
    base_group = is_qa_group ? GROUP[1:(end - 3)] : GROUP

    if contains(base_group, "OrdinaryDiffEq") || base_group == "ImplicitDiscreteSolve" || base_group == "SimpleImplicitDiscreteSolve"
        Pkg.activate(joinpath(dirname(@__DIR__), "lib", base_group))
        # Set QA_ONLY env var to tell sublibrary tests whether to run only QA tests
        withenv("ODEDIFFEQ_TEST_GROUP" => (is_qa_group ? "QA" : "FUNCTIONAL")) do
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
        @time @safetestset "Differentiation Trait Tests" include("interface/differentiation_traits_tests.jl")
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
        @time @safetestset "Linear Nonlinear Solver Tests" include("interface/linear_nonlinear_tests.jl")
        @time @safetestset "Linear Solver Tests" include("interface/linear_solver_test.jl")
        @time @safetestset "Linear Solver Split ODE Tests" include("interface/linear_solver_split_ode_test.jl")
        @time @safetestset "AutoSparse Detection Tests" include("interface/autosparse_detection_tests.jl")
        @time @safetestset "Enum Tests" include("interface/enums.jl")
        @time @safetestset "CheckInit Tests" include("interface/checkinit_tests.jl")
        @time @safetestset "Get du Tests" include("interface/get_du.jl")
        @time @safetestset "Mass Matrix Tests" include("interface/mass_matrix_tests.jl")
        @time @safetestset "W-Operator prototype tests" include("interface/wprototype_tests.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceIII" || GROUP == "Interface")
        @time @safetestset "Derivative Utilities Tests" include("interface/utility_tests.jl")
        @time @safetestset "stats Tests" include("interface/stats_tests.jl")
        @time @safetestset "No Index Tests" include("interface/noindex_tests.jl")
        @time @safetestset "Events + DAE addsteps Tests" include("interface/event_dae_addsteps.jl")
        @time @safetestset "No Jac Tests" include("interface/nojac.jl")
        @time @safetestset "Units Tests" include("interface/units_tests.jl")
        @time @safetestset "Non-Full Diagonal Sparsity Tests" include("interface/nonfulldiagonal_sparse.jl")
        @time @safetestset "DEVerbosity Tests" include("interface/verbosity.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceIV" || GROUP == "Interface")
        @time @safetestset "Autodiff Error Tests" include("interface/autodiff_error_tests.jl")
        @time @safetestset "Ambiguity Tests" include("interface/ambiguity_tests.jl")
        @time @safetestset "Precision Mixing Tests" include("interface/precision_mixing.jl")
        @time @safetestset "Sized Matrix Tests" include("interface/sized_matrix_tests.jl")
        @time @safetestset "Second Order with First Order Solver Tests" include("interface/second_order_with_first_order_solvers.jl")
    end

    if !is_APPVEYOR && (GROUP == "All" || GROUP == "InterfaceV" || GROUP == "Interface")
        @time @safetestset "Interpolation Derivative Error Tests" include("interface/interpolation_derivative_error_tests.jl")
        @time @safetestset "GPU AutoDiff Interface Tests" include("interface/gpu_autodiff_interface_tests.jl")
        @time @safetestset "DAE Initialization Tests" include("interface/dae_initialization_tests.jl")
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

    # Don't run ModelingToolkit tests on prerelease
    if !is_APPVEYOR && GROUP == "ModelingToolkit" && isempty(VERSION.prerelease)
        activate_modelingtoolkit_env()
        @time @safetestset "NLStep Tests" include("modelingtoolkit/nlstep_tests.jl")
        @time @safetestset "Jacobian Tests" include("modelingtoolkit/jacobian_tests.jl")
        @time @safetestset "Preconditioner Tests" include("modelingtoolkit/preconditioners.jl")
        @time @safetestset "DAE Initialize Integration" include("modelingtoolkit/dae_initialize_integration.jl")
    end

    if !is_APPVEYOR && GROUP == "Downstream"
        activate_downstream_env()
        @time @safetestset "DelayDiffEq Tests" include("downstream/delaydiffeq.jl")
        @time @safetestset "Measurements Tests" include("downstream/measurements.jl")
        @time @safetestset "Sparse Diff Tests" include("downstream/sparsediff_tests.jl")
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

    if !is_APPVEYOR && GROUP == "Multithreading"
        @time @safetestset "Extrapolation Tests" include("multithreading/ode_extrapolation_tests.jl")
    end

    if !is_APPVEYOR && GROUP == "GPU"
        activate_gpu_env()
        @time @safetestset "Simple GPU" begin
            import OrdinaryDiffEqCore
            include(
                joinpath(
                    dirname(pathof(OrdinaryDiffEqCore.DiffEqBase)), "..",
                    "test/gpu/simple_gpu.jl"
                )
            )
        end
        @time @safetestset "Autoswitch GPU" include("gpu/autoswitch.jl")
        @time @safetestset "Linear LSRK GPU" include("gpu/linear_lsrk.jl")
        @time @safetestset "Linear Exponential GPU" include("gpu/linear_exp.jl")
        @time @safetestset "Reaction-Diffusion Stiff Solver GPU" include("gpu/reaction_diffusion_stiff.jl")
        @time @safetestset "Scalar indexing bug bypass" include("gpu/hermite_test.jl")
        @time @safetestset "RKIP Semilinear PDE GPU" include("gpu/rkip_semilinear_pde.jl")
        @time @safetestset "simple dae on GPU" include("gpu/simple_dae.jl")
        @time @safetestset "BDF solvers GPU" include("gpu/bdf_solvers.jl")
    end

    if !is_APPVEYOR && GROUP == "QA"
        @time @safetestset "Quality Assurance Tests" include("qa/qa_tests.jl")
    end
end # @time
