using SafeTestsets
using Pkg

function activate_gpu_env()
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

function activate_nopre_env()
    Pkg.activate("nopre")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    return Pkg.instantiate()
end

const LONGER_TESTS = false

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV, "APPVEYOR")

@time begin
    if TEST_GROUP == "ALL" || TEST_GROUP == "Core" || TEST_GROUP == "Interface1"
        @time @safetestset "Default Solver Tests" begin
            include("default_solver_test.jl")
        end
        @time @safetestset "First Rand Tests" begin
            include("first_rand_test.jl")
        end
        @time @safetestset "Inference Tests" begin
            include("inference_test.jl")
        end
        @time @safetestset "Linear RODE Tests" begin
            include("rode_linear_tests.jl")
        end
        @time @safetestset "Complex Number Tests" begin
            include("complex_tests.jl")
        end
        @time @safetestset "Static Array Tests" begin
            include("static_array_tests.jl")
        end
        @time @safetestset "StaticArrays Callback Tests (issue #1129)" begin
            include("staticarrays_callback_test.jl")
        end
        @time @safetestset "Noise Type Tests" begin
            include("noise_type_test.jl")
        end
        @time @safetestset "Mass matrix tests" begin
            include("mass_matrix_tests.jl")
        end
        @time @safetestset "Outofplace Arrays Tests" begin
            include("outofplace_arrays.jl")
        end
        @time @safetestset "tdir Tests" begin
            include("tdir_tests.jl")
        end
        @time @safetestset "tstops Tests" begin
            include("tstops_tests.jl")
        end
        @time @safetestset "Callable tstops Tests" begin
            include("callable_tstops_tests.jl")
        end
        @time @safetestset "Integrator RNG Tests" begin
            include("rng_integrator_tests.jl")
        end
        @time @safetestset "JumpProblem kwarg forwarding Tests" begin
            include("jump_kwarg_forwarding_tests.jl")
        end
        @time @safetestset "saveat Tests" begin
            include("saveat_tests.jl")
        end
        @time @safetestset "Oval2" begin
            include("oval2_test.jl")
        end
        @time @safetestset "Solver Reversal Tests" begin
            include("reversal_tests.jl")
        end
    end

    if TEST_GROUP == "ALL" || TEST_GROUP == "Interface2"
        @time @safetestset "Basic Tau Leaping Tests" begin
            include("tau_leaping.jl")
        end
        @time @safetestset "Linear SDE Tests" begin
            include("sde/sde_linear_tests.jl")
        end
        @time @safetestset "Two-dimensional Linear SDE Tests" begin
            include("sde/sde_twodimlinear_tests.jl")
        end
        @time @safetestset "Element-wise Tolerances Tests" begin
            include("tolerances_tests.jl")
        end
        @time @safetestset "Zeroed Noise Tests" begin
            include("zeroed_noise_test.jl")
        end
        @time @safetestset "Scalar Tests" begin
            include("scalar_noise.jl")
        end
        @time @safetestset "Stiffness Detection Test" begin
            include("stiffness_detection_test.jl")
        end
        @time @safetestset "Adaptive SDE Linear Tests" begin
            include("adaptive/sde_linearadaptive_tests.jl")
        end
        @time @safetestset "Inplace RODESolution Interpolation Tests" begin
            include("inplace_interpolation.jl")
        end
    end

    if TEST_GROUP == "ALL" || TEST_GROUP == "Interface3"
        @time @safetestset "Composite Tests" begin
            include("composite_algorithm_test.jl")
        end
        @time @safetestset "Events Tests" begin
            include("events_test.jl")
        end
        @time @safetestset "Cache Tests" begin
            include("cache_test.jl")
        end
        @time @safetestset "Adaptive Complex Mean Test" begin
            include("adaptive/sde_complex_adaptive_mean_test.jl")
        end
        @time @safetestset "Utility Tests" begin
            include("utility_tests.jl")
        end
        @time @safetestset "Implicit Time Parameter Tests" begin
            include("implicit_time_parameter_test.jl")
        end
        @time @safetestset "Non-diagonal SDE Tests" begin
            include("nondiagonal_tests.jl")
        end
        @time @safetestset "Non-diagonal EulerHeun sparse alloc" begin
            include("nondiag_noise_eulerheun_test.jl")
        end
        @time @safetestset "No Index Tests" begin
            include("noindex_tests.jl")
        end
        @time @safetestset "Multiple Dimension Linear Adaptive Test" begin
            include("adaptive/sde_twodimlinearadaptive_tests.jl")
        end
        @time @safetestset "Autostepsize Test" begin
            include("adaptive/sde_autostepsize_test.jl")
        end
        @time @safetestset "Additive Lorenz Attractor Test" begin
            include("adaptive/sde_lorenzattractor_tests.jl")
        end
        @time @safetestset "Stochastic iterated integrals" begin
            include("levy_areas.jl")
        end
        @time @safetestset "Sparse Diff Tests" begin
            include("sparsediff_tests.jl")
        end
    end

    if !is_APPVEYOR && (TEST_GROUP == "ALL" || TEST_GROUP == "AlgConvergence")
        @time @safetestset "Convergence Tests" begin
            include("sde/sde_convergence_tests.jl")
        end
        @time @safetestset "Dynamical SDE Tests" begin
            include("sde/sde_dynamical.jl")
        end
    end

    if !is_APPVEYOR && (TEST_GROUP == "ALL" || TEST_GROUP == "NoncommutativeConvergence")
        @time @safetestset "Noncommutative Noise Tests" begin
            include("noncommutative_tests.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "AlgConvergence2"
        @time @safetestset "IIF Convergence Tests" begin
            include("iif_methods.jl")
        end
        @time @safetestset "Commutative Noise Methods Tests" begin
            include("commutative_tests.jl")
        end
        @time @safetestset "Multivariate Geometric Tests" begin
            include("multivariate_geometric.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "AlgConvergence3"
        @time @safetestset "Rossler Order Tests" begin
            include("sde/sde_rosslerorder_tests.jl")
        end
        @time @safetestset "ODE Convergence Regression Tests" begin
            include("ode_convergence_regression.jl")
        end
        @time @safetestset "Additive SDE Tests" begin
            include("sde/sde_additive_tests.jl")
        end
        @time @safetestset "Stratonovich Convergence Tests" begin
            include("stratonovich_convergence_tests.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "WeakConvergence1"
        @time @safetestset "Multidimensional IIP Weak Convergence Tests" begin
            include("weak_convergence/multidim_iip_weak.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "OOPWeakConvergence"
        @time @safetestset "OOP Weak Convergence Tests" begin
            include("weak_convergence/oop_weak.jl")
        end
        @time @safetestset "Additive Weak Convergence Tests" begin
            include("weak_convergence/additive_weak.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "IIPWeakConvergence"
        @time @safetestset "IIP Weak Convergence Tests" begin
            include("weak_convergence/iip_weak.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "WeakAdaptiveGPU"
        activate_gpu_env()
        @time @safetestset "GPU Weak adaptive step size scalar noise SDE" begin
            include("gpu/sde_weak_scalar_adaptive_gpu.jl")
        end
        @time @safetestset "GPU Weak adaptive" begin
            include("gpu/sde_weak_adaptive_gpu.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "Multithreaded"
        @time @safetestset "Mulithreaded Jump Thread Safety Tests" begin
            include("multithreaded_jump_test.jl")
        end
    end

    if TEST_GROUP == "ALL" || TEST_GROUP == "Allocations"
        @time @safetestset "Allocation Tests" begin
            include("alloc_tests.jl")
        end
    end

    if !is_APPVEYOR && TEST_GROUP == "StaticAnalysis"
        activate_nopre_env()
        @time @safetestset "JET Static Analysis Tests" begin
            include("nopre/jet_tests.jl")
        end
    end
end
