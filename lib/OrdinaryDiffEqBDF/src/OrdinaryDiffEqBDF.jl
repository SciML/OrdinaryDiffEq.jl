module OrdinaryDiffEqBDF

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, @unpack, unwrap_alg,
                           calculate_residuals, alg_extrapolates,
                           OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqNewtonAlgorithm,
                           AbstractController, DEFAULT_PRECS,
                           CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache,
                           isfsal, full_cache,
                           constvalue, isadaptive, error_constant,
                           has_special_newton_error,
                           trivial_limiter!,
                           issplit, qsteady_min_default, qsteady_max_default,
                           get_current_alg_order, get_current_adaptive_order,
                           default_controller, stepsize_controller!,
                           step_accept_controller!,
                           step_reject_controller!, post_newton_controller!,
                           u_modified!, DAEAlgorithm, _unwrap_val, DummyController,
                           get_fsalfirstlast, generic_solver_docstring
using OrdinaryDiffEqSDIRK: ImplicitEulerConstantCache, ImplicitEulerCache

using TruncatedStacktraces, MuladdMacro, MacroTools, FastBroadcast, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA
using LinearAlgebra: mul!, I
using ArrayInterface
import OrdinaryDiffEqCore
using OrdinaryDiffEqDifferentiation: UJacobianWrapper
using OrdinaryDiffEqNonlinearSolve: NLNewton, du_alias_or_new, build_nlsolver,
                                    nlsolve!, nlsolvefail, isnewton, markfirststage!,
                                    set_new_W!, DIRK, compute_step!, COEFFICIENT_MULTISTEP

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("bdf_utils.jl")
include("bdf_caches.jl")
include("dae_caches.jl")
include("controllers.jl")
include("dae_perform_step.jl")
include("bdf_perform_step.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [FBDF()]
    prob_list = []

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]))
    end

    for prob in prob_list, solver in solver_list
        solve(prob, solver)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export ABDF2, QNDF1, QBDF1, QNDF2, QBDF2, QNDF, QBDF, FBDF,
       SBDF, SBDF2, SBDF3, SBDF4, MEBDF2, IMEXEuler, IMEXEulerARK,
       DABDF2, DImplicitEuler, DFBDF

end
