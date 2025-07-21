"""
$(DocStringExtensions.README)
"""
module OrdinaryDiffEq

import Reexport
using Reexport: @reexport
@reexport using DiffEqBase

# Explicit imports for functions that are re-exported
import DiffEqBase: DiffEqBase, init, solve, solve!, step!, addsteps!, savevalues!, terminate!

import OrdinaryDiffEqCore: OrdinaryDiffEqCore,
                           CompositeAlgorithm,
                           ShampineCollocationInit, BrownFullBasicInit, NoInit,
                           copyat_or_push!,
                           du_cache, full_cache, isfsal, ode_interpolant, u_cache,
                           AutoSwitch,
                           @unpack,
                           @cache

export CompositeAlgorithm, ShampineCollocationInit, BrownFullBasicInit, NoInit,
       AutoSwitch

import OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqDifferentiation: OrdinaryDiffEqDifferentiation

using OrdinaryDiffEqNonlinearSolve: NLNewton, NLAnderson, NLFunctional,
                                    NonlinearSolveAlg
export NLNewton, NLAnderson, NLFunctional, NonlinearSolveAlg

using OrdinaryDiffEqExtrapolation: AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
                                   ImplicitEulerExtrapolation, ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
                                   ImplicitEulerBarycentricExtrapolation
export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
       ImplicitEulerExtrapolation,
       ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
       ImplicitEulerBarycentricExtrapolation

using OrdinaryDiffEqStabilizedRK: ROCK2, ROCK4, RKC, ESERK4, ESERK5, SERK2
export ROCK2, ROCK4, RKC, ESERK4, ESERK5, SERK2

using OrdinaryDiffEqStabilizedIRK: IRKC
export IRKC

using OrdinaryDiffEqLowStorageRK: ORK256, CarpenterKennedy2N54, SHLDDRK64, HSLDDRK64, DGLDDRK73_C, DGLDDRK84_C,
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

using OrdinaryDiffEqSSPRK: SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
                          SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H,
                          SSPRK33, KYKSSPRK42, KYK2014DGSSPRK_3S2
export SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
       SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H,
       SSPRK33, KYKSSPRK42, KYK2014DGSSPRK_3S2

using OrdinaryDiffEqFeagin: Feagin10, Feagin12, Feagin14
export Feagin10, Feagin12, Feagin14

using OrdinaryDiffEqSymplecticRK: SymplecticEuler, VelocityVerlet, VerletLeapfrog, LeapfrogDriftKickDrift,
                                   PseudoVerletLeapfrog, McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
                                   CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10
export SymplecticEuler, VelocityVerlet, VerletLeapfrog, LeapfrogDriftKickDrift,
       PseudoVerletLeapfrog, McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
       CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

using OrdinaryDiffEqRKN: Nystrom4, FineRKN4, FineRKN5, Nystrom4VelocityIndependent,
                         Nystrom5VelocityIndependent,
                         IRKN3, IRKN4, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, DPRKN12, ERKN4, ERKN5, ERKN7,
                         RKN4
export Nystrom4, FineRKN4, FineRKN5, Nystrom4VelocityIndependent,
       Nystrom5VelocityIndependent,
       IRKN3, IRKN4, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, DPRKN12, ERKN4, ERKN5, ERKN7,
       RKN4

using OrdinaryDiffEqVerner: Vern6, Vern7, Vern8, Vern9
export Vern6, Vern7, Vern8, Vern9

using OrdinaryDiffEqSDIRK: ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22,
                           Kvaerno3, KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4,
                           Kvaerno5, KenCarp4, KenCarp47, KenCarp5, KenCarp58, ESDIRK54I8L2SA, SFSDIRK4,
                           SFSDIRK5, CFNLIRK3, SFSDIRK6, SFSDIRK7, SFSDIRK8,
                           ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2, ESDIRK659L2SA
export ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22,
       Kvaerno3, KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4,
       Kvaerno5, KenCarp4, KenCarp47, KenCarp5, KenCarp58, ESDIRK54I8L2SA, SFSDIRK4,
       SFSDIRK5, CFNLIRK3, SFSDIRK6, SFSDIRK7, SFSDIRK8,
       ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2, ESDIRK659L2SA

using OrdinaryDiffEqBDF: ABDF2, QNDF1, QBDF1, QNDF2, QBDF2, QNDF, QBDF, FBDF,
                         SBDF2, SBDF3, SBDF4, MEBDF2, IMEXEuler, IMEXEulerARK,
                         DImplicitEuler, DABDF2, DFBDF
export ABDF2, QNDF1, QBDF1, QNDF2, QBDF2, QNDF, QBDF, FBDF,
       SBDF2, SBDF3, SBDF4, MEBDF2, IMEXEuler, IMEXEulerARK,
       DImplicitEuler, DABDF2, DFBDF

using OrdinaryDiffEqTsit5: Tsit5, AutoTsit5
export Tsit5, AutoTsit5

using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
                                Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
                                Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
                                RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
                                ROS3PRL2, ROK4a,
                                ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7
export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7

using OrdinaryDiffEqDefault: DefaultODEAlgorithm
export DefaultODEAlgorithm

using OrdinaryDiffEqFIRK: RadauIIA3, RadauIIA5, RadauIIA9
export RadauIIA3, RadauIIA5, RadauIIA9

using OrdinaryDiffEqQPRK: QPRK98
export QPRK98

using OrdinaryDiffEqPDIRK: PDIRK44
export PDIRK44

using OrdinaryDiffEqPRK: KuttaPRK2p5
export KuttaPRK2p5


using OrdinaryDiffEqLowOrderRK: Euler, SplitEuler, Heun, Ralston, Midpoint, RK4,
                                BS3, OwrenZen3, OwrenZen4, OwrenZen5, BS5,
                                DP5, Anas5, RKO65, FRK65, RKM, MSRK5, MSRK6,
                                PSRK4p7q6, PSRK3p5q4, PSRK3p6q5, Stepanov5, SIR54,
                                Alshina2, Alshina3, Alshina6, AutoDP5
export Euler, SplitEuler, Heun, Ralston, Midpoint, RK4,
       BS3, OwrenZen3, OwrenZen4, OwrenZen5, BS5,
       DP5, Anas5, RKO65, FRK65, RKM, MSRK5, MSRK6,
       PSRK4p7q6, PSRK3p5q4, PSRK3p6q5, Stepanov5, SIR54,
       Alshina2, Alshina3, Alshina6, AutoDP5

using OrdinaryDiffEqFunctionMap: FunctionMap
export FunctionMap

using OrdinaryDiffEqAdamsBashforthMoulton: AB3, AB4, AB5, ABM32, ABM43, ABM54, VCAB3,
                                          VCAB4, VCAB5, VCABM3, VCABM4, VCABM5, VCABM
export AB3, AB4, AB5, ABM32, ABM43, ABM54, VCAB3,
       VCAB4, VCAB5, VCABM3, VCABM4, VCABM5, VCABM

using OrdinaryDiffEqNordsieck: AN5, JVODE, JVODE_Adams, JVODE_BDF
export AN5, JVODE, JVODE_Adams, JVODE_BDF

using OrdinaryDiffEqExplicitRK: ExplicitRK, constructDormandPrince
export ExplicitRK

using OrdinaryDiffEqLinear: MagnusMidpoint, LinearExponential, MagnusLeapfrog, LieEuler, CayleyEuler,
                           MagnusGauss4, MagnusNC6, MagnusGL6, MagnusGL8, MagnusNC8, MagnusGL4,
                           MagnusAdapt4, RKMK2, RKMK4, LieRK4, CG2, CG3, CG4a
export MagnusMidpoint, LinearExponential, MagnusLeapfrog, LieEuler, CayleyEuler,
       MagnusGauss4, MagnusNC6, MagnusGL6, MagnusGL8, MagnusNC8, MagnusGL4,
       MagnusAdapt4, RKMK2, RKMK4, LieRK4, CG2, CG3, CG4a

using OrdinaryDiffEqIMEXMultistep: CNAB2, CNLF2
export CNAB2, CNLF2

using OrdinaryDiffEqExponentialRK: LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
                                   EPIRK4s3B,
                                   EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43
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
