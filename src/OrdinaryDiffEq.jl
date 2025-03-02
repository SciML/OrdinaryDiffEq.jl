"""
$(DocStringExtensions.README)
"""
module OrdinaryDiffEq

using Reexport
@reexport using DiffEqBase

import OrdinaryDiffEqCore: OrdinaryDiffEqCore,
                           trivial_limiter!, CompositeAlgorithm, alg_order,
                           ShampineCollocationInit, BrownFullBasicInit, NoInit,
                           set_new_W!, set_W_γdt!, get_W, isfirstcall, isfirststage,
                           isJcurrent, get_new_W_γdt_cutoff,
                           DIRK, COEFFICIENT_MULTISTEP, NORDSIECK_MULTISTEP, GLM,
                           MethodType, Divergence, VerySlowConvergence,
                           SlowConvergence, Convergence, FastConvergence, NLStatus,
                           TryAgain, AbstractNLSolverCache,
                           AbstractNLSolverAlgorithm, AbstractNLSolver,
                           handle_discontinuities!, copyat_or_push!,
                           du_cache, full_cache, isfsal, ode_interpolant, u_cache,
                           AutoSwitch, has_discontinuity,
                           first_discontinuity, pop_discontinuity!, _vec, loopfooter!,
                           _reshape, perform_step!,
                           _ode_addsteps!, get_current_alg_autodiff, default_controller,
                           isstandard,
                           ispredictive, beta2_default, beta1_default, gamma_default,
                           qmin_default,
                           qmax_default, qsteady_min_default, qsteady_max_default,
                           stepsize_controller!,
                           accept_step_controller, step_accept_controller!,
                           step_reject_controller!,
                           DummyController, issplit, calculate_residuals,
                           calculate_residuals!,
                           nlsolve_f, unwrap_cache, ode_addsteps!, get_chunksize,
                           handle_callback_modifiers!,
                           unwrap_alg, apply_step!, initialize_tstops, uses_uprev,
                           initialize_saveat,
                           isimplicit, initialize_d_discontinuities, isdtchangeable,
                           _searchsortedfirst,
                           _searchsortedlast,
                           @unpack, ismultistep, DEFAULT_PRECS, isautoswitch,
                           get_chunksize_int,
                           _unwrap_val, alg_autodiff, concrete_jac, alg_difftype,
                           standardtag,
                           alg_extrapolates, alg_maximum_order, alg_adaptive_order,
                           OrdinaryDiffEqCompositeAlgorithm, initialize_callbacks!,
                           PredictiveController,
                           get_differential_vars, alg_cache, AutoSwitchCache,
                           InterpolationData,
                           DEOptions, OrdinaryDiffEqAlgorithm, @cache, fsal_typeof,
                           OrdinaryDiffEqCache,
                           OrdinaryDiffEqAdaptiveAlgorithm, handle_dt!,
                           ode_determine_initdt,
                           loopheader!, OrdinaryDiffEqConstantCache, _loopfooter!,
                           isadaptive,
                           OrdinaryDiffEqMutableCache, current_interpolant!,
                           is_mass_matrix_alg,
                           False, True, _savevalues!, postamble!, recursivefill!,
                           _change_t_via_interpolation!, ODEIntegrator, _ode_interpolant!,
                           current_interpolant, resize_nlsolver!, _ode_interpolant,
                           handle_tstop!, _postamble!, update_uprev!, resize_J_W!,
                           DAEAlgorithm, get_fsalfirstlast, strip_cache,
                           Sequential, BaseThreads, PolyesterThreads

export CompositeAlgorithm, ShampineCollocationInit, BrownFullBasicInit, NoInit
AutoSwitch

import OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqDifferentiation: _alg_autodiff, resize_grad_config!, dolinsolve,
                                     wrapprecs, UJacobianWrapper, build_jac_config,
                                     WOperator, FirstAutodiffJacError, calc_J!, calc_W!,
                                     calc_J, calc_W, jacobian2W!, isnewton

using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLAnderson, NLFunctional, nlsolvefail,
                                    initial_η, NonlinearSolveAlg, compute_step!, NLSolver,
                                    nlsolve!, resize_jac_config!, anderson!, build_nlsolver,
                                    markfirststage!, anderson
export NLNewton, NLAnderson, NLFunctional, NonlinearSolveAlg

using OrdinaryDiffEqExtrapolation
export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
       ImplicitEulerExtrapolation,
       ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
       ImplicitEulerBarycentricExtrapolation

using OrdinaryDiffEqStabilizedRK
export ROCK2, ROCK4, RKC, ESERK4, ESERK5, SERK2

using OrdinaryDiffEqStabilizedIRK
export IRKC

using OrdinaryDiffEqLowStorageRK
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
       RK46NL, SHLDDRK_2N, SHLDDRK52

using OrdinaryDiffEqSSPRK
export SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
       SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H,
       SSPRK33, KYKSSPRK42, KYK2014DGSSPRK_3S2

using OrdinaryDiffEqFeagin
export Feagin10, Feagin12, Feagin14

using OrdinaryDiffEqSymplecticRK
export SymplecticEuler, VelocityVerlet, VerletLeapfrog, LeapfrogDriftKickDrift,
       PseudoVerletLeapfrog, McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
       CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

using OrdinaryDiffEqRKN
export Nystrom4, FineRKN4, FineRKN5, Nystrom4VelocityIndependent,
       Nystrom5VelocityIndependent,
       IRKN3, IRKN4, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, DPRKN12, ERKN4, ERKN5, ERKN7,
       RKN4

using OrdinaryDiffEqVerner
export Vern6, Vern7, Vern8, Vern9

using OrdinaryDiffEqSDIRK
import OrdinaryDiffEqSDIRK: ImplicitEulerConstantCache, ImplicitEulerCache
export ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22,
       Kvaerno3, KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4,
       Kvaerno5, KenCarp4, KenCarp47, KenCarp5, KenCarp58, ESDIRK54I8L2SA, SFSDIRK4,
       SFSDIRK5, CFNLIRK3, SFSDIRK6, SFSDIRK7, SFSDIRK8, Kvaerno5, KenCarp4, KenCarp5,
       SFSDIRK4, SFSDIRK5, CFNLIRK3, SFSDIRK6,
       SFSDIRK7, SFSDIRK8, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2, ESDIRK659L2SA

using OrdinaryDiffEqBDF
export ABDF2, QNDF1, QBDF1, QNDF2, QBDF2, QNDF, QBDF, FBDF,
       SBDF2, SBDF3, SBDF4, MEBDF2, IMEXEuler, IMEXEulerARK,
       DImplicitEuler, DABDF2, DFBDF

using OrdinaryDiffEqTsit5
export Tsit5, AutoTsit5
import OrdinaryDiffEqTsit5: Tsit5ConstantCache, Tsit5Cache

using OrdinaryDiffEqRosenbrock
export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7
import OrdinaryDiffEqRosenbrock: RosenbrockMutableCache

using OrdinaryDiffEqDefault
export DefaultODEAlgorithm

using OrdinaryDiffEqFIRK
export RadauIIA3, RadauIIA5, RadauIIA9

using OrdinaryDiffEqQPRK
export QPRK98

using OrdinaryDiffEqPDIRK
export PDIRK44

using OrdinaryDiffEqPRK
export KuttaPRK2p5

using OrdinaryDiffEqHighOrderRK
export TanYam7, DP8, PFRK87, TsitPap8
using OrdinaryDiffEqHighOrderRK: DP8ConstantCache, DP8Cache

using OrdinaryDiffEqLowOrderRK
export Euler, SplitEuler, Heun, Ralston, Midpoint, RK4,
       BS3, OwrenZen3, OwrenZen4, OwrenZen5, BS5,
       DP5, Anas5, RKO65, FRK65, RKM, MSRK5, MSRK6,
       PSRK4p7q6, PSRK3p5q4, PSRK3p6q5, Stepanov5, SIR54,
       Alshina2, Alshina3, Alshina6, AutoDP5
using OrdinaryDiffEqLowOrderRK: BS3Cache, BS3ConstantCache, RK4ConstantCache, RK4Cache

using OrdinaryDiffEqFunctionMap
export FunctionMap

using OrdinaryDiffEqAdamsBashforthMoulton
export AB3, AB4, AB5, ABM32, ABM43, ABM54, VCAB3,
       VCAB4, VCAB5, VCABM3, VCABM4, VCABM5, VCABM

using OrdinaryDiffEqNordsieck
export AN5, JVODE, JVODE_Adams, JVODE_BDF

using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructDormandPrince
export ExplicitRK

using OrdinaryDiffEqLinear
export MagnusMidpoint, LinearExponential, MagnusLeapfrog, LieEuler, CayleyEuler,
       MagnusGauss4, MagnusNC6, MagnusGL6, MagnusGL8, MagnusNC8, MagnusGL4,
       MagnusAdapt4, RKMK2, RKMK4, LieRK4, CG2, CG3, CG4a

using OrdinaryDiffEqIMEXMultistep
export CNAB2, CNLF2

using OrdinaryDiffEqExponentialRK
export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
       EPIRK4s3B,
       EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43

import PrecompileTools
import Preferences
import DocStringExtensions

const DEPRECATED_ADDSTEPS = true

#General Functions
export solve, solve!, init, step!

#Callback Necessary
export addsteps!, ode_interpolant, terminate!, savevalues!, copyat_or_push!, isfsal

export constructDormandPrince

# Reexport the Alg Types

export CompositeAlgorithm

export AutoSwitch, AutoTsit5, AutoDP5,
       AutoVern6, AutoVern7, AutoVern8, AutoVern9

import OrdinaryDiffEqCore: IController, PIController, PIDController
export IController, PIController, PIDController
end # module
