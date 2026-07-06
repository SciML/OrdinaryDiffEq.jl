module OrdinaryDiffEqSDIRK

# These OrdinaryDiffEqCore functions are extended with SDIRK-specific methods, so
# they must be brought in with `import` (not `using`) to allow method definitions.
import OrdinaryDiffEqCore: perform_step!,
    alg_extrapolates,
    alg_cache, full_cache,
    isesdirk, issplit,
    ssp_coefficient, get_fsalfirstlast
# OrdinaryDiffEqCore names used (called/referenced) but not extended here.
using OrdinaryDiffEqCore: unwrap_alg,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    OrdinaryDiffEqNewtonAlgorithm,
    CompiledFloats,
    constvalue,
    trivial_limiter!,
    generic_solver_docstring,
    _fixup_ad, current_extrapolant!, Predictor,
    isnewton, get_W, set_new_W!, COEFFICIENT_MULTISTEP,
    find_algebraic_vars_eqs
export Predictor
using TruncatedStacktraces: @truncate_stacktrace
using MuladdMacro: MuladdMacro, @muladd
using MacroTools: MacroTools
using FastBroadcast: FastBroadcast, @..
using RecursiveArrayTools: RecursiveArrayTools, recursivefill!
# `alg_order` is owned by SciMLBase and extended here, so it needs `import`.
import SciMLBase: alg_order
using SciMLBase: SciMLBase, SplitFunction, ODEProblem, _vec, _reshape, _unwrap_val
# `initialize!` is owned by DiffEqBase and extended here, so it needs `import`;
# `calculate_residuals`/`calculate_residuals!` are only called.
import DiffEqBase: initialize!
using DiffEqBase: calculate_residuals, calculate_residuals!
using LinearAlgebra: mul!, I
import OrdinaryDiffEqCore

using OrdinaryDiffEqDifferentiation: dolinsolve
using OrdinaryDiffEqNonlinearSolve: du_alias_or_new, markfirststage!, build_nlsolver,
    nlsolve!, nlsolvefail,
    NLNewton
import ADTypes: AutoForwardDiff
using CommonSolve: solve

using Reexport: Reexport, @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("sdirk_caches.jl")
include("kencarp_kvaerno_caches.jl")
include("sdirk_perform_step.jl")
include("kencarp_kvaerno_perform_step.jl")
include("sdirk_tableaus.jl")
include("imex_tableaus.jl")
include("generic_imex_perform_step.jl")

export ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22,
    Kvaerno3, KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4,
    Kvaerno5, KenCarp4, KenCarp47, KenCarp5, KenCarp58, ESDIRK54I8L2SA, SFSDIRK4,
    SFSDIRK5, CFNLIRK3, SFSDIRK6, SFSDIRK7, SFSDIRK8, Kvaerno5, KenCarp4, KenCarp5,
    SFSDIRK4, SFSDIRK5, CFNLIRK3, SFSDIRK6,
    SFSDIRK7, SFSDIRK8, ESDIRK325L2SA, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2,
    ESDIRK659L2SA,
    ARS343, ARS222, ARS232, ARS443,
    IMEXSSP222, IMEXSSP2322, IMEXSSP3332, IMEXSSP3433, BHR553

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [TRBDF2(), KenCarp4()]
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

# Cross-sublibrary IMEX cache/tableau API that other OrdinaryDiffEq solver
# sublibraries (e.g. OrdinaryDiffEqBDF) reference to build IMEX methods on top
# of the ESDIRK-IMEX machinery. Marked public so those references are
# recognized as a supported extension API rather than internal access.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Expr(
            :public,
            :ESDIRKIMEXCache, :ESDIRKIMEXConstantCache,
            :ImplicitEulerESDIRKIMEXTableau
        )
    )
end

end
