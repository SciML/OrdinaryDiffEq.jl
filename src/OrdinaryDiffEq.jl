"""
$(DocStringExtensions.README)
"""
module OrdinaryDiffEq

if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

using DocStringExtensions
using Reexport
@reexport using DiffEqBase

using Logging

using MuladdMacro, SparseArrays, FastClosures

using LinearAlgebra

using PrecompileTools

import FillArrays: Trues, Falses

# Interfaces
import DiffEqBase: solve!, step!, initialize!, isadaptive

# Internal utils
import DiffEqBase: ODE_DEFAULT_NORM,
                   ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE,
                   ODE_DEFAULT_UNSTABLE_CHECK

import SciMLOperators: SciMLOperators, AbstractSciMLOperator, AbstractSciMLScalarOperator,
                       MatrixOperator, FunctionOperator,
                       update_coefficients, update_coefficients!, DEFAULT_UPDATE_FUNC,
                       isconstant

using DiffEqBase: DEIntegrator

import RecursiveArrayTools: chain, recursivecopy!

using SimpleUnPack, RecursiveArrayTools,
      DataStructures, ArrayInterface, ArrayInterface

import TruncatedStacktraces

using ExponentialUtilities

# Required by temporary fix in not in-place methods with 12+ broadcasts
# `MVector` is used by Nordsieck forms
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

# Integrator Interface
import DiffEqBase: resize!, deleteat!, addat!, full_cache, user_cache, u_cache, du_cache,
                   resize_non_user_cache!, deleteat_non_user_cache!, addat_non_user_cache!,
                   terminate!, get_du, get_dt, get_proposed_dt, set_proposed_dt!,
                   u_modified!, savevalues!,
                   add_tstop!, has_tstop, first_tstop, pop_tstop!,
                   add_saveat!, set_reltol!,
                   set_abstol!, postamble!, last_step_failed,
                   isautodifferentiable

using DiffEqBase: check_error!, @def, _vec, _reshape

using FastBroadcast: @.., True, False

using SciMLBase: NoInit, _unwrap_val

import SciMLBase: alg_order

import DiffEqBase: calculate_residuals,
                   calculate_residuals!, unwrap_cache,
                   @tight_loop_macros,
                   islinear, timedepentdtmin

import Polyester
using MacroTools, Adapt
import ADTypes: AutoFiniteDiff, AutoForwardDiff

using SciMLStructures: canonicalize, Tunable, isscimlstructure

const CompiledFloats = Union{Float32, Float64}
import Preferences

abstract type AbstractNLSolverCache end
abstract type AbstractNLSolverAlgorithm end
abstract type AbstractNLSolver{algType, iip} end

function set_new_W! end
function set_W_γdt! end
function get_W end

DEFAULT_PRECS(W, du, u, p, t, newW, Plprev, Prprev, solverdata) = nothing, nothing
isdiscretecache(cache) = false

include("doc_utils.jl")
include("misc_utils.jl")

include("algorithms.jl")
include("composite_algs.jl")

include("alg_utils.jl")

include("caches/basic_caches.jl")

include("integrators/type.jl")
include("integrators/controllers.jl")
include("integrators/integrator_interface.jl")
include("integrators/integrator_utils.jl")
include("cache_utils.jl")
include("initialize_dae.jl")

include("perform_step/composite_perform_step.jl")

include("dense/generic_dense.jl")

include("iterator_interface.jl")
include("solve.jl")
include("initdt.jl")
include("interp_func.jl")

include("../lib/OrdinaryDiffEqDifferentiation/src/OrdinaryDiffEqDifferentiation.jl")
import ..OrdinaryDiffEqDifferentiation

include("../lib/OrdinaryDiffEqNonlinearSolve/src/OrdinaryDiffEqNonlinearSolve.jl")
using ..OrdinaryDiffEqNonlinearSolve

include("../lib/OrdinaryDiffEqExtrapolation/src/OrdinaryDiffEqExtrapolation.jl")
using ..OrdinaryDiffEqExtrapolation
export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
       ImplicitEulerExtrapolation,
       ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
       ImplicitEulerBarycentricExtrapolation

include("../lib/OrdinaryDiffEqStabilizedRK/src/OrdinaryDiffEqStabilizedRK.jl")
using ..OrdinaryDiffEqStabilizedRK
export ROCK2, ROCK4, RKC, ESERK4, ESERK5, SERK2

include("../lib/OrdinaryDiffEqStabilizedIRK/src/OrdinaryDiffEqStabilizedIRK.jl")
using ..OrdinaryDiffEqStabilizedIRK
export IRKC

include("../lib/OrdinaryDiffEqLowStorageRK/src/OrdinaryDiffEqLowStorageRK.jl")
using ..OrdinaryDiffEqLowStorageRK
export ORK256, CarpenterKennedy2N54, SHLDDRK64, HSLDDRK64, DGLDDRK73_C, DGLDDRK84_C,
       DGLDDRK84_F, NDBLSRK124, NDBLSRK134, NDBLSRK144,
       CFRLDDRK64, TSLDDRK74, CKLLSRK43_2, CKLLSRK54_3C, CKLLSRK95_4S, CKLLSRK95_4C,
       CKLLSRK95_4M,
       CKLLSRK54_3C_3R, CKLLSRK54_3M_3R, CKLLSRK54_3N_3R, CKLLSRK85_4C_3R, CKLLSRK85_4M_3R,
       CKLLSRK85_4P_3R,
       CKLLSRK54_3N_4R, CKLLSRK54_3M_4R, CKLLSRK65_4M_4R, CKLLSRK85_4FM_4R, CKLLSRK75_4M_5R,
       ParsaniKetchesonDeconinck3S32, ParsaniKetchesonDeconinck3S82,
       ParsaniKetchesonDeconinck3S53, ParsaniKetchesonDeconinck3S173,
       ParsaniKetchesonDeconinck3S94, ParsaniKetchesonDeconinck3S184,
       ParsaniKetchesonDeconinck3S105, ParsaniKetchesonDeconinck3S205,
       RDPK3Sp35, RDPK3SpFSAL35, RDPK3Sp49, RDPK3SpFSAL49, RDPK3Sp510, RDPK3SpFSAL510,
       KYK2014DGSSPRK_3S2, RK46NL

include("../lib/OrdinaryDiffEqSSPRK/src/OrdinaryDiffEqSSPRK.jl")
using ..OrdinaryDiffEqSSPRK
export SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
       SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H,
       SSPRK33, SHLDDRK_2N, KYKSSPRK42, SHLDDRK52

include("../lib/OrdinaryDiffEqFeagin/src/OrdinaryDiffEqFeagin.jl")
using ..OrdinaryDiffEqFeagin
export Feagin10, Feagin12, Feagin14

include("../lib/OrdinaryDiffEqSymplecticRK/src/OrdinaryDiffEqSymplecticRK.jl")
using ..OrdinaryDiffEqSymplecticRK
export SymplecticEuler, VelocityVerlet, VerletLeapfrog, PseudoVerletLeapfrog,
       McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
       CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

include("../lib/OrdinaryDiffEqRKN/src/OrdinaryDiffEqRKN.jl")
using ..OrdinaryDiffEqRKN
export Nystrom4, FineRKN4, FineRKN5, Nystrom4VelocityIndependent,
       Nystrom5VelocityIndependent,
       IRKN3, IRKN4, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, DPRKN12, ERKN4, ERKN5, ERKN7,
       RKN4

include("../lib/OrdinaryDiffEqVerner/src/OrdinaryDiffEqVerner.jl")
using ..OrdinaryDiffEqVerner
export Vern6, Vern7, Vern8, Vern9

include("../lib/OrdinaryDiffEqSDIRK/src/OrdinaryDiffEqSDIRK.jl")
using ..OrdinaryDiffEqSDIRK
import .OrdinaryDiffEqSDIRK: ImplicitEulerConstantCache, ImplicitEulerCache
export ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22,
       Kvaerno3, KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4,
       Kvaerno5, KenCarp4, KenCarp47, KenCarp5, KenCarp58, ESDIRK54I8L2SA, SFSDIRK4,
       SFSDIRK5, CFNLIRK3, SFSDIRK6, SFSDIRK7, SFSDIRK8, Kvaerno5, KenCarp4, KenCarp5,
       SFSDIRK4, SFSDIRK5, CFNLIRK3, SFSDIRK6,
       SFSDIRK7, SFSDIRK8, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2, ESDIRK659L2SA

include("../lib/OrdinaryDiffEqBDF/src/OrdinaryDiffEqBDF.jl")
using ..OrdinaryDiffEqBDF
export ABDF2, QNDF1, QBDF1, QNDF2, QBDF2, QNDF, QBDF, FBDF,
       SBDF2, SBDF3, SBDF4, MEBDF2, IMEXEuler, IMEXEulerARK,
       DImplicitEuler, DABDF2, DFBDF

include("../lib/OrdinaryDiffEqTsit5/src/OrdinaryDiffEqTsit5.jl")
using ..OrdinaryDiffEqTsit5
export Tsit5, AutoTsit5
import .OrdinaryDiffEqTsit5: Tsit5ConstantCache, Tsit5Cache

include("../lib/OrdinaryDiffEqRosenbrock/src/OrdinaryDiffEqRosenbrock.jl")
using ..OrdinaryDiffEqRosenbrock
export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7
import ..OrdinaryDiffEqRosenbrock: RosenbrockMutableCache

include("../lib/OrdinaryDiffEqDefault/src/OrdinaryDiffEqDefault.jl")
using ..OrdinaryDiffEqDefault
export DefaultODEAlgorithm

include("../lib/OrdinaryDiffEqFIRK/src/OrdinaryDiffEqFIRK.jl")
using ..OrdinaryDiffEqFIRK
export RadauIIA3, RadauIIA5, RadauIIA9

include("../lib/OrdinaryDiffEqQPRK/src/OrdinaryDiffEqQPRK.jl")
using ..OrdinaryDiffEqQPRK
export QPRK98

include("../lib/OrdinaryDiffEqPDIRK/src/OrdinaryDiffEqPDIRK.jl")
using ..OrdinaryDiffEqPDIRK
export PDIRK44

include("../lib/OrdinaryDiffEqPRK/src/OrdinaryDiffEqPRK.jl")
using ..OrdinaryDiffEqPRK
export KuttaPRK2p5

include("../lib/OrdinaryDiffEqHighOrderRK/src/OrdinaryDiffEqHighOrderRK.jl")
using ..OrdinaryDiffEqHighOrderRK
export TanYam7, DP8, PFRK87, TsitPap8
using ..OrdinaryDiffEqHighOrderRK: DP8ConstantCache, DP8Cache

include("../lib/OrdinaryDiffEqLowOrderRK/src/OrdinaryDiffEqLowOrderRK.jl")
using ..OrdinaryDiffEqLowOrderRK
export Euler, SplitEuler, Heun, Ralston, Midpoint, RK4,
       BS3, OwrenZen3, OwrenZen4, OwrenZen5, BS5,
       DP5, Anas5, RKO65, FRK65, RKM, MSRK5, MSRK6,
       PSRK4p7q6, PSRK3p5q4, PSRK3p6q5, Stepanov5, SIR54,
       Alshina2, Alshina3, Alshina6, AutoDP5
using ..OrdinaryDiffEqLowOrderRK: BS3Cache, BS3ConstantCache, RK4ConstantCache, RK4Cache

include("../lib/OrdinaryDiffEqFunctionMap/src/OrdinaryDiffEqFunctionMap.jl")
using ..OrdinaryDiffEqFunctionMap
export FunctionMap

include("../lib/OrdinaryDiffEqAdamsBashforthMoulton/src/OrdinaryDiffEqAdamsBashforthMoulton.jl")
using ..OrdinaryDiffEqAdamsBashforthMoulton
export AB3, AB4, AB5, ABM32, ABM43, ABM54, VCAB3,
       VCAB4, VCAB5, VCABM3, VCABM4, VCABM5, VCABM

include("../lib/OrdinaryDiffEqNordsieck/src/OrdinaryDiffEqNordsieck.jl")
using ..OrdinaryDiffEqNordsieck
export AN5, JVODE, JVODE_Adams, JVODE_BDF

include("../lib/OrdinaryDiffEqExplicitRK/src/OrdinaryDiffEqExplicitRK.jl")
using ..OrdinaryDiffEqExplicitRK
using ..OrdinaryDiffEqExplicitRK: constructDormandPrince
export ExplicitRK

include("../lib/OrdinaryDiffEqLinear/src/OrdinaryDiffEqLinear.jl")
using ..OrdinaryDiffEqLinear
export MagnusMidpoint, LinearExponential, MagnusLeapfrog, LieEuler, CayleyEuler,
       MagnusGauss4, MagnusNC6, MagnusGL6, MagnusGL8, MagnusNC8, MagnusGL4,
       MagnusAdapt4, RKMK2, RKMK4, LieRK4, CG2, CG3, CG4a

include("../lib/OrdinaryDiffEqIMEXMultistep/src/OrdinaryDiffEqIMEXMultistep.jl")
using ..OrdinaryDiffEqIMEXMultistep
export CNAB2, CNLF2

include("../lib/OrdinaryDiffEqExponentialRK/src/OrdinaryDiffEqExponentialRK.jl")
using ..OrdinaryDiffEqExponentialRK
export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
       EPIRK4s3B,
       EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43

PrecompileTools.@compile_workload begin
    function lorenz(du, u, p, t)
        du[1] = 10.0(u[2] - u[1])
        du[2] = u[1] * (28.0 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
    end

    function lorenz_oop(u, p, t)
        [10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
    end

    nonstiff = [
        Tsit5(), Vern7()
    ]

    stiff = [Rosenbrock23(),
        Rodas5P(),
        FBDF()
    ]

    default_ode = [
        DefaultODEAlgorithm(autodiff = false)
    ]

    default_autodiff_ode = [
        DefaultODEAlgorithm()
    ]

    autoswitch = [
        AutoTsit5(Rosenbrock23(autodiff = false)),
        AutoTsit5(TRBDF2(autodiff = false)),
        AutoVern7(Rodas5P(autodiff = false)),
        AutoVern7(KenCarp47(autodiff = false))
    ]

    low_storage = [
        SSPRK43(), RDPK3SpFSAL35(), RDPK3SpFSAL49()
    ]

    low_storage_nonadaptive = [
        CarpenterKennedy2N54(williamson_condition = false)
    ]

    solver_list = []
    solver_list_nonadaptive = []

    if Preferences.@load_preference("PrecompileDefault", true)
        append!(solver_list, default_ode)
    end

    if Preferences.@load_preference("PrecompileAutodiffDefault", true)
        append!(solver_list, default_autodiff_ode)
    end

    if Preferences.@load_preference("PrecompileNonStiff", false)
        append!(solver_list, nonstiff)
    end

    if Preferences.@load_preference("PrecompileStiff", false)
        append!(solver_list, stiff)
    end

    if Preferences.@load_preference("PrecompileAutoSwitch", false)
        append!(solver_list, autoswitch)
    end

    if Preferences.@load_preference("PrecompileLowStorage", false)
        append!(solver_list, low_storage)
        append!(solver_list_nonadaptive, low_storage_nonadaptive)
    end

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

    for prob in prob_list, solver in solver_list_nonadaptive
        solve(prob, solver; dt = 0.5)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

const DEPRECATED_ADDSTEPS = true

#General Functions
export solve, solve!, init, step!

#Callback Necessary
export addsteps!, ode_interpolant, terminate!, savevalues!, copyat_or_push!, isfsal

export constructDormandPrince

# Reexport the Alg Types

export CompositeAlgorithm

export SHLDDRK52, SHLDDRK_2N

export AutoSwitch, AutoTsit5, AutoDP5,
       AutoVern6, AutoVern7, AutoVern8, AutoVern9

export ShampineCollocationInit, BrownFullBasicInit, NoInit

export NLNewton, NLAnderson, NLFunctional, NonlinearSolveAlg

export IController, PIController, PIDController
end # module
