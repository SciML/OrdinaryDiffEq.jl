module OrdinaryDiffEqBDF

import OrdinaryDiffEqCore: perform_step!, unwrap_alg,
    alg_extrapolates,
    default_controller, IController,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    OrdinaryDiffEqNewtonAlgorithm,
    AbstractController,
    alg_cache, @cache,
    isfsal, full_cache,
    constvalue, error_constant,
    has_special_newton_error,
    trivial_limiter!,
    issplit, qmax_default, gamma_default,
    qsteady_min_default, qsteady_max_default,
    get_current_alg_order, get_current_adaptive_order,
    stepsize_controller!,
    step_accept_controller!,
    step_reject_controller!, post_newton_controller!,
    get_EEst,
    setup_controller_cache, get_qmax, get_gamma, get_qsteady_min, get_qsteady_max,
    get_failfactor, CommonControllerOptions, resolve_basic, _resolved_QT,
    AbstractControllerCache,
    DAEAlgorithm,
    get_fsalfirstlast, generic_solver_docstring, _fixup_ad,
    _ode_interpolant, _ode_interpolant!, has_stiff_interpolation,
    _ode_addsteps!, DerivativeOrderNotPossibleError, set_discontinuity,
    DIRK, COEFFICIENT_MULTISTEP, isnewton, set_new_W!,
    find_algebraic_vars_eqs
import SciMLBase: alg_order, isadaptive, _unwrap_val
import DiffEqBase: calculate_residuals, calculate_residuals!, initialize!
using OrdinaryDiffEqSDIRK: ESDIRKIMEXConstantCache, ESDIRKIMEXCache,
    ImplicitEulerESDIRKIMEXTableau

using TruncatedStacktraces: @truncate_stacktrace
using MuladdMacro: @muladd
using FastBroadcast: @..
using RecursiveArrayTools: recursivefill!
using LinearAlgebra: mul!, I, Diagonal
import ArrayInterface
import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: default_controller

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

using OrdinaryDiffEqNonlinearSolve: NLNewton, du_alias_or_new, build_nlsolver,
    nlsolve!, nlsolvefail, markfirststage!
import ADTypes: AutoForwardDiff

using Reexport: Reexport, @reexport
import SciMLBase
using SciMLBase: ODEProblem, ODEFunction, derivative_discontinuity!, solve
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

# Index-1 mass matrix DAE (Robertson). The pure-ODE problems in the workload below
# never reach the algebraic variable detection or the DAE initialization solve, so
# without this every mass matrix FBDF user pays that inference at first solve.
function precompile_mm_dae(du, u, p, t)
    du[1] = -0.04u[1] + 1.0e4 * u[2] * u[3]
    du[2] = 0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3]
    return du[3] = u[1] + u[2] + u[3] - 1.0
end

PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [FBDF()]
    prob_list = []

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileMassMatrixDAE", true)
        mm_dae = ODEFunction(
            precompile_mm_dae, mass_matrix = Diagonal([1.0, 1.0, 0.0])
        )
        push!(prob_list, ODEProblem(mm_dae, [1.0, 0.0, 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(mm_dae, [1.0, 0.0, 0.0], (0.0, 1.0), Float64[]))
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

    if Preferences.@load_preference("PrecompileAutoDePSpecialize", false)
        push!(
            prob_list,
            ODEProblem{true, OrdinaryDiffEqCore.AutoDePSpecialize}(
                OrdinaryDiffEqCore.lorenz_p, [1.0; 0.0; 0.0],
                (0.0, 1.0), OrdinaryDiffEqCore.lorenz_p_params
            )
        )
        push!(
            prob_list,
            ODEProblem{true, OrdinaryDiffEqCore.AutoDePSpecialize}(
                OrdinaryDiffEqCore.lorenz_pref, [1.0; 0.0; 0.0],
                (0.0, 1.0), OrdinaryDiffEqCore.lorenz_pref_params
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
