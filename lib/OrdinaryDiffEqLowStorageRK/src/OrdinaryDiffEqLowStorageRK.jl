module OrdinaryDiffEqLowStorageRK

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order, calculate_residuals!,
    beta2_default, beta1_default, gamma_default,
    initialize!, perform_step!, unwrap_alg,
    calculate_residuals, ssp_coefficient,
    OrdinaryDiffEqAlgorithm, ispredictive,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqAdaptiveAlgorithm, uses_uprev,
    default_controller_v7,
    legacy_default_controller, NewPIDController, PIDController,
    alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
    constvalue, _unwrap_val,
    trivial_limiter!, perform_step!, initialize!,
    explicit_rk_docstring, get_fsalfirstlast
using FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools, Adapt
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA
import Static: False
import RecursiveArrayTools: recursive_unitless_bottom_eltype
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

include("arrayfuse.jl")
include("algorithms.jl")
include("alg_utils.jl")
include("low_storage_rk_caches.jl")
include("low_storage_rk_perform_step.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = []
    solver_list_nonadaptive = []
    prob_list = []

    low_storage = [
        RDPK3SpFSAL35(), RDPK3SpFSAL49(),
    ]

    low_storage_nonadaptive = [
        CarpenterKennedy2N54(williamson_condition = false),
    ]

    if Preferences.@load_preference("PrecompileLowStorage", false)
        append!(solver_list, low_storage)
        append!(solver_list_nonadaptive, low_storage_nonadaptive)
    end

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

    for prob in prob_list, solver in solver_list_nonadaptive
        solve(prob, solver; dt = 0.5)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export ORK256, CarpenterKennedy2N54, SHLDDRK64, HSLDDRK64, DGLDDRK73_C, DGLDDRK84_C,
    DGLDDRK84_F, NDBLSRK124, NDBLSRK134, NDBLSRK144,
    CFRLDDRK64, TSLDDRK74, CKLLSRK43_2, CKLLSRK54_3C,
    CKLLSRK95_4S, CKLLSRK95_4C, CKLLSRK95_4M,
    CKLLSRK54_3C_3R, CKLLSRK54_3M_3R, CKLLSRK54_3N_3R, CKLLSRK85_4C_3R, CKLLSRK85_4M_3R,
    CKLLSRK85_4P_3R,
    CKLLSRK54_3N_4R, CKLLSRK54_3M_4R, CKLLSRK65_4M_4R, CKLLSRK85_4FM_4R, CKLLSRK75_4M_5R,
    ParsaniKetchesonDeconinck3S32, ParsaniKetchesonDeconinck3S82,
    ParsaniKetchesonDeconinck3S53, ParsaniKetchesonDeconinck3S173,
    ParsaniKetchesonDeconinck3S94, ParsaniKetchesonDeconinck3S184,
    ParsaniKetchesonDeconinck3S105, ParsaniKetchesonDeconinck3S205,
    RDPK3Sp35, RDPK3SpFSAL35, RDPK3Sp49, RDPK3SpFSAL49, RDPK3Sp510, RDPK3SpFSAL510,
    RK46NL, SHLDDRK_2N, SHLDDRK52
end
