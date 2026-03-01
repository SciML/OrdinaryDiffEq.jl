module OrdinaryDiffEqBDF

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
    initialize!, perform_step!, unwrap_alg,
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
    stepsize_controller!,
    step_accept_controller!,
    step_reject_controller!, post_newton_controller!,
    u_modified!, DAEAlgorithm, _unwrap_val, DummyController,
    get_fsalfirstlast, generic_solver_docstring, _bool_to_ADType,
    _process_AD_choice,
    _ode_interpolant, _ode_interpolant!, has_stiff_interpolation,
    _ode_addsteps!, DerivativeOrderNotPossibleError,
    initdt_alg
using OrdinaryDiffEqSDIRK: ImplicitEulerConstantCache, ImplicitEulerCache

using TruncatedStacktraces: @truncate_stacktrace
using MuladdMacro: @muladd
using MacroTools: @capture
using FastBroadcast: @..
using RecursiveArrayTools: recursivefill!
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA
using LinearAlgebra: mul!, I
import ArrayInterface
using ArrayInterface: ismutable
import OrdinaryDiffEqCore

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        import OrdinaryDiffEqCore: default_controller_v7,
            legacy_default_controller
    end
else
    @eval begin
        import OrdinaryDiffEqCore: default_controller
    end
end

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.10"
    @eval begin
        import OrdinaryDiffEqCore: get_current_qmax
    end
else
    @eval begin
        # Fallback for older OrdinaryDiffEqCore: no first-step qmax behavior
        @inline get_current_qmax(integrator, qmax) = qmax
    end
end

using OrdinaryDiffEqDifferentiation: UJacobianWrapper
using OrdinaryDiffEqNonlinearSolve: NLNewton, du_alias_or_new, build_nlsolver,
    nlsolve!, nlsolvefail, isnewton, markfirststage!,
    set_new_W!, DIRK, compute_step!, COEFFICIENT_MULTISTEP,
    NonlinearSolveAlg
import ADTypes: AutoForwardDiff, AutoFiniteDiff, AbstractADType

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("bdf_utils.jl")
include("stald.jl")
include("bdf_caches.jl")
include("dae_caches.jl")
include("controllers.jl")
include("dae_perform_step.jl")
include("bdf_perform_step.jl")
include("interp_func.jl")
include("bdf_interpolants.jl")
include("stiff_addsteps.jl")

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
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]
            )
        )
    end

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)
            )
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]
            )
        )
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0))
        )
        push!(
            prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(
                lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]
            )
        )
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
