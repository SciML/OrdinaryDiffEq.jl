abstract type OrdinaryDiffEqAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

abstract type OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end
abstract type OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end

abstract type OrdinaryDiffEqImplicitAlgorithm{CS,AD} <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqNewtonAlgorithm{CS,AD} <:  OrdinaryDiffEqImplicitAlgorithm{CS,AD} end
abstract type OrdinaryDiffEqRosenbrockAlgorithm{CS,AD} <:  OrdinaryDiffEqImplicitAlgorithm{CS,AD} end
const NewtonAlgorithm = Union{OrdinaryDiffEqNewtonAlgorithm,OrdinaryDiffEqNewtonAdaptiveAlgorithm}
const RosenbrockAlgorithm = Union{OrdinaryDiffEqRosenbrockAlgorithm,OrdinaryDiffEqRosenbrockAdaptiveAlgorithm}

abstract type OrdinaryDiffEqExponentialAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqAdaptiveExponentialAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
const ExponentialAlgorithm = Union{OrdinaryDiffEqExponentialAlgorithm,OrdinaryDiffEqAdaptiveExponentialAlgorithm}

abstract type OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end

struct FunctionMap{scale_by_time} <: OrdinaryDiffEqAlgorithm end
FunctionMap(;scale_by_time=false) = FunctionMap{scale_by_time}()

###############################################################################

# RK methods

struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
  tableau::TabType
end
ExplicitRK(;tableau=ODE_DEFAULT_TABLEAU) = ExplicitRK(tableau)

@inline trivial_limiter!(u, f, t) = nothing

struct Euler <: OrdinaryDiffEqAlgorithm end
struct KuttaPRK2p5 <: OrdinaryDiffEqAlgorithm
  threading::Bool
end
KuttaPRK2p5(;threading=true) = KuttaPRK2p5(threading)

struct AitkenNeville <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
  max_order::Int
  min_order::Int
  init_order::Int
  threading::Bool
end
AitkenNeville(;max_order=10,min_order=1,init_order=5,threading=true) = AitkenNeville(max_order,min_order,init_order,threading)

struct ImplicitEulerExtrapolation{CS,AD,F,FDT} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  max_order::Int
  min_order::Int
  init_order::Int
  threading::Bool
  diff_type::FDT
end

function ImplicitEulerExtrapolation(;chunk_size=0,autodiff=true,
    diff_type=Val{:forward},linsolve=DEFAULT_LINSOLVE,
    max_order=10,min_order=1,init_order=5,threading=true)
    if threading
      @warn "Threading in `ImplicitEulerExtrapolation` is currently disabled. Thus `threading` has been changed from `true` to `false`."
      threading = false
    end
    ImplicitEulerExtrapolation{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(
      linsolve,max_order,min_order,init_order,threading,diff_type)
end

struct ExtrapolationMidpointDeuflhard <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  threading::Bool
end
function ExtrapolationMidpointDeuflhard(;min_order=1,init_order=5, max_order=10, sequence = :harmonic, threading = true)
  # Enforce 1 <=  min_order <= init_order <= max_order:
  n_min = max(1,min_order)
  n_init = max(n_min,init_order)
  n_max = max(n_init,max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ExtrapolationMidpointDeuflhard` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
  end

  # Warn user if sequence has been changed:
  if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
    @warn "The `sequence` given to the `ExtrapolationMidpointDeuflhard` algorithm
       is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
       Thus it has been changed
      :$(sequence) --> :harmonic"
    sequence = :harmonic
  end

  # Initialize algorithm
  ExtrapolationMidpointDeuflhard(n_min,n_init,n_max,sequence,threading)
end

struct ImplicitDeuflhardExtrapolation{CS,AD,F,FDT} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  diff_type::FDT
end
function ImplicitDeuflhardExtrapolation(;chunk_size=0,autodiff=true,
  linsolve=DEFAULT_LINSOLVE,diff_type=Val{:forward},
  min_order=1,init_order=5,max_order=10,sequence = :harmonic)
  # Enforce 1 <=  min_order <= init_order <= max_order:
  n_min = max(1,min_order)
  n_init = max(n_min,init_order)
  n_max = max(n_init,max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ImplicitDeuflhardExtrapolation` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
  end

  # Warn user if sequence has been changed:
  if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
    @warn "The `sequence` given to the `ImplicitDeuflhardExtrapolation` algorithm
       is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
       Thus it has been changed
      :$(sequence) --> :harmonic"
    sequence = :harmonic
  end

  # Initialize algorithm
  ImplicitDeuflhardExtrapolation{chunk_size, autodiff,
      typeof(linsolve), typeof(diff_type)}(linsolve,n_min,n_init,n_max,sequence,diff_type)
end

struct ExtrapolationMidpointHairerWanner <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  threading::Bool
end
function ExtrapolationMidpointHairerWanner(;min_order=2,init_order=5, max_order=10, sequence = :harmonic, threading = true)
  # Enforce 2 <=  min_order
  # and min_order + 1 <= init_order <= max_order - 1:
  n_min = max(2, min_order)
  n_init = max(n_min + 1, init_order)
  n_max = max(n_init + 1, max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ExtrapolationMidpointHairerWanner` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
  end

  # Warn user if sequence has been changed:
  if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
    @warn "The `sequence` given to the `ExtrapolationMidpointHairerWanner` algorithm
       is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
       Thus it has been changed
      :$(sequence) --> :harmonic"
    sequence = :harmonic
  end

  # Initialize algorithm
  ExtrapolationMidpointHairerWanner(n_min,n_init,n_max,sequence,threading)
end

struct ImplicitHairerWannerExtrapolation{CS,AD,F,FDT} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  diff_type::FDT
end
function ImplicitHairerWannerExtrapolation(;chunk_size=0,autodiff=true,
  linsolve=DEFAULT_LINSOLVE,diff_type=Val{:forward},
  min_order=2,init_order=5,max_order=10,sequence = :harmonic)
  # Enforce 2 <=  min_order
  # and min_order + 1 <= init_order <= max_order - 1:
  n_min = max(2, min_order)
  n_init = max(n_min + 1, init_order)
  n_max = max(n_init + 1, max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ImplicitHairerWannerExtrapolation` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
  end

  # Warn user if sequence has been changed:
  if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
    @warn "The `sequence` given to the `ImplicitHairerWannerExtrapolation` algorithm
       is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
       Thus it has been changed
      :$(sequence) --> :harmonic"
    sequence = :harmonic
  end

  # Initialize algorithm
  ImplicitHairerWannerExtrapolation{chunk_size, autodiff,
      typeof(linsolve), typeof(diff_type)}(linsolve,n_min,n_init,n_max,
      sequence,diff_type)
end

struct RK46NL <: OrdinaryDiffEqAlgorithm end
struct Heun <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Ralston <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Midpoint <: OrdinaryDiffEqAdaptiveAlgorithm end
struct RK4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct RKM <: OrdinaryDiffEqAlgorithm end
struct Anas5{T} <: OrdinaryDiffEqAlgorithm
  w::T
end
Anas5(; w=1) = Anas5(w)

struct ORK256 <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  ORK256(;williamson_condition=true) = new(williamson_condition)
end
struct CarpenterKennedy2N54 <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  CarpenterKennedy2N54(;williamson_condition=true) = new(williamson_condition)
end
struct SHLDDRK52 <: OrdinaryDiffEqAlgorithm end
struct SHLDDRK_2N <: OrdinaryDiffEqAlgorithm end
struct HSLDDRK64 <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  HSLDDRK64(;williamson_condition=true) = new(williamson_condition)
end
struct DGLDDRK73_C <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  DGLDDRK73_C(;williamson_condition=true) = new(williamson_condition)
end
struct DGLDDRK84_C <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  DGLDDRK84_C(;williamson_condition=true) = new(williamson_condition)
end
struct DGLDDRK84_F <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  DGLDDRK84_F(;williamson_condition=true) = new(williamson_condition)
end
struct NDBLSRK124 <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  NDBLSRK124(;williamson_condition=true) = new(williamson_condition)
end
struct NDBLSRK134 <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  NDBLSRK134(;williamson_condition=true) = new(williamson_condition)
end
struct NDBLSRK144 <: OrdinaryDiffEqAlgorithm
  williamson_condition::Bool
  NDBLSRK144(;williamson_condition=true) = new(williamson_condition)
end
struct CFRLDDRK64 <: OrdinaryDiffEqAlgorithm end
struct TSLDDRK74 <: OrdinaryDiffEqAlgorithm end
struct CKLLSRK43_2 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3C <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK95_4S <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK95_4C <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK95_4M <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3C_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3M_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3N_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4C_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4M_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4P_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3N_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3M_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK65_4M_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4FM_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK75_4M_5R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct ParsaniKetchesonDeconinck3S32 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S82 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S53 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S173 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S94 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S184 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S105 <: OrdinaryDiffEqAlgorithm end
struct ParsaniKetchesonDeconinck3S205 <: OrdinaryDiffEqAlgorithm end
struct KYK2014DGSSPRK_3S2 <: OrdinaryDiffEqAlgorithm end

struct RKO65 <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

struct SSPRK22{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK22(stage_limiter! = trivial_limiter!) = SSPRK22(stage_limiter!, trivial_limiter!)
struct SSPRK33{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK33(stage_limiter! = trivial_limiter!) = SSPRK33(stage_limiter!, trivial_limiter!)
struct SSPRK53{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
struct KYKSSPRK42 <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
SSPRK53(stage_limiter! = trivial_limiter!) = SSPRK53(stage_limiter!, trivial_limiter!)
struct SSPRK53_2N1{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK53_2N1(stage_limiter! = trivial_limiter!) = SSPRK53_2N1(stage_limiter!, trivial_limiter!)
struct SSPRK53_2N2{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK53_2N2(stage_limiter! = trivial_limiter!) = SSPRK53_2N2(stage_limiter!, trivial_limiter!)
struct SSPRK53_H{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK53_H(stage_limiter! = trivial_limiter!) = SSPRK53_H(stage_limiter!, trivial_limiter!)
struct SSPRK63{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK63(stage_limiter! = trivial_limiter!) = SSPRK63(stage_limiter!, trivial_limiter!)
struct SSPRK73{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK73(stage_limiter! = trivial_limiter!) = SSPRK73(stage_limiter!, trivial_limiter!)
struct SSPRK83{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK83(stage_limiter! = trivial_limiter!) = SSPRK83(stage_limiter!, trivial_limiter!)
struct SSPRK432{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK432(stage_limiter! = trivial_limiter!) = SSPRK432(stage_limiter!, trivial_limiter!)
struct SSPRKMSVS43{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRKMSVS43(stage_limiter! = trivial_limiter!) = SSPRKMSVS43(stage_limiter!, trivial_limiter!)
struct SSPRKMSVS32{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRKMSVS32(stage_limiter! = trivial_limiter!) = SSPRKMSVS32(stage_limiter!, trivial_limiter!)
struct SSPRK932{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK932(stage_limiter! = trivial_limiter!) = SSPRK932(stage_limiter!, trivial_limiter!)
struct SSPRK54{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK54(stage_limiter! = trivial_limiter!) = SSPRK54(stage_limiter!, trivial_limiter!)
struct SSPRK104{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK104(stage_limiter! = trivial_limiter!) = SSPRK104(stage_limiter!, trivial_limiter!)

struct OwrenZen3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct OwrenZen4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct OwrenZen5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct BS3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Tsit5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP8 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct TanYam7 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct TsitPap8 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Feagin10 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Feagin12 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Feagin14 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct BS5 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  BS5(;lazy=true) = new(lazy)
end
struct Vern6 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern6(;lazy=true) = new(lazy)
end
struct Vern7 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern7(;lazy=true) = new(lazy)
end
struct Vern8 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern8(;lazy=true) = new(lazy)
end
struct Vern9 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern9(;lazy=true) = new(lazy)
end
struct FRK65{T} <: OrdinaryDiffEqAdaptiveAlgorithm
  omega::T
  FRK65(omega=0.0) = new{typeof(omega)}(omega)
end

struct PFRK87{T} <: OrdinaryDiffEqAdaptiveAlgorithm
  omega::T
  PFRK87(omega=0.0) = new{typeof(omega)}(omega)
end


################################################################################

# Symplectic methods

struct SymplecticEuler <: OrdinaryDiffEqAlgorithm end
struct VelocityVerlet <: OrdinaryDiffEqAlgorithm end
struct VerletLeapfrog <: OrdinaryDiffEqAlgorithm end
struct PseudoVerletLeapfrog <: OrdinaryDiffEqAlgorithm end
struct McAte2 <: OrdinaryDiffEqAlgorithm end
struct Ruth3 <: OrdinaryDiffEqAlgorithm end
struct McAte3 <: OrdinaryDiffEqAlgorithm end
struct CandyRoz4 <: OrdinaryDiffEqAlgorithm end
struct McAte4 <: OrdinaryDiffEqAlgorithm end
struct CalvoSanz4 <: OrdinaryDiffEqAlgorithm end
struct McAte42 <: OrdinaryDiffEqAlgorithm end
struct McAte5 <: OrdinaryDiffEqAlgorithm end
struct Yoshida6 <: OrdinaryDiffEqAlgorithm end
struct KahanLi6 <: OrdinaryDiffEqAlgorithm end
struct McAte8 <: OrdinaryDiffEqAlgorithm end
struct KahanLi8 <: OrdinaryDiffEqAlgorithm end
struct SofSpa10 <: OrdinaryDiffEqAlgorithm end

# Nyström methods

struct IRKN3 <: OrdinaryDiffEqAlgorithm end
struct Nystrom4 <: OrdinaryDiffEqAlgorithm end
struct Nystrom4VelocityIndependent <: OrdinaryDiffEqAlgorithm end
struct IRKN4 <: OrdinaryDiffEqAlgorithm end
struct Nystrom5VelocityIndependent <: OrdinaryDiffEqAlgorithm end
struct DPRKN6 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DPRKN8 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DPRKN12 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct ERKN4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct ERKN5 <: OrdinaryDiffEqAdaptiveAlgorithm end

################################################################################

# Adams Bashforth and Adams moulton methods

struct AB3 <: OrdinaryDiffEqAlgorithm end
struct AB4 <: OrdinaryDiffEqAlgorithm end
struct AB5 <: OrdinaryDiffEqAlgorithm end
struct ABM32 <: OrdinaryDiffEqAlgorithm end
struct ABM43 <: OrdinaryDiffEqAlgorithm end
struct ABM54 <: OrdinaryDiffEqAlgorithm end


# Variable Step Size Adams methods

struct VCAB3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct VCAB4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct VCAB5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct VCABM3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct VCABM4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct VCABM5 <: OrdinaryDiffEqAdaptiveAlgorithm end

# Variable Order and Variable Step Size Adams methods

struct VCABM <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm end

# IMEX Multistep methods

struct CNAB2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CNAB2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CNAB2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct CNLF2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CNLF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CNLF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct QNDF1{CS,AD,F,F2,FDT,κType} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  kappa::κType
  controller::Symbol
end
QNDF1(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                  extrapolant=:linear,kappa = -0.1850,
                 controller = :Standard) =
                 QNDF1{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                 typeof(kappa)}(
                 linsolve,nlsolve,diff_type,extrapolant,kappa,controller)

QBDF1(;kwargs...) = QNDF1(;kappa=0,kwargs...)

struct QNDF2{CS,AD,F,F2,FDT,κType} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  kappa::κType
  controller::Symbol
end
QNDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                  extrapolant=:linear,kappa = -1//9,
                 controller = :Standard) =
                 QNDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                 typeof(kappa)}(
                 linsolve,nlsolve,diff_type,extrapolant,kappa,controller)

QBDF2(;kwargs...) = QNDF2(;kappa=0,kwargs...)


struct QNDF{CS,AD,F,F2,FDT,K,T,κType} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  kappa::κType
  controller::Symbol
end
Base.@pure QNDF(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),κ=nothing,tol=nothing,
                extrapolant=:linear,kappa=promote(-0.1850,-1//9,-0.0823,-0.0415,0),
                controller = :Standard) =
                QNDF{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                typeof(κ),typeof(tol),typeof(kappa)}(
                linsolve,nlsolve,diff_type,κ,tol,extrapolant,kappa,controller)

Base.@pure QBDF(;kwargs...) = QNDF(;kappa=tuple(0//1,0//1,0//1,0//1,0//1),kwargs...)

struct SBDF{CS,AD,F,F2,FDT,K,T} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  order::Int
end
SBDF(order;chunk_size=0,autodiff=true,diff_type=Val{:forward},
     linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),κ=nothing,tol=nothing,
     extrapolant=:linear) =
     SBDF{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
     typeof(κ),typeof(tol)}(
     linsolve,nlsolve,diff_type,κ,tol,extrapolant,order)
IMEXEuler(;kwargs...) = SBDF(1;kwargs...)
SBDF2(;kwargs...) = SBDF(2;kwargs...)
SBDF3(;kwargs...) = SBDF(3;kwargs...)
SBDF4(;kwargs...) = SBDF(4;kwargs...)

# Adams/BDF methods in Nordsieck forms
struct AN5   <: OrdinaryDiffEqAdaptiveAlgorithm end
struct JVODE{bType,aType} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
  algorithm::Symbol
  bias1::bType
  bias2::bType
  bias3::bType
  addon::aType
end

JVODE(algorithm=:Adams;bias1=6, bias2=6,bias3=10,
                 addon=1//10^6) = JVODE(algorithm,bias1,bias2,bias3,addon)
JVODE_Adams(;kwargs...) = JVODE(:Adams;kwargs...)
JVODE_BDF(;kwargs...) = JVODE(:BDF;kwargs...)

# ROCK methods
struct ROCK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
  min_stages::Int
  max_stages::Int
  eigen_est::E
end
ROCK2(;min_stages=0,max_stages=200,eigen_est=nothing) = ROCK2(min_stages,max_stages,eigen_est)

struct ROCK4{E} <: OrdinaryDiffEqAdaptiveAlgorithm
  min_stages::Int
  max_stages::Int
  eigen_est::E
end
ROCK4(;min_stages=0,max_stages=152,eigen_est=nothing) = ROCK4(min_stages,max_stages,eigen_est)

# SERK methods
for Alg in [:ESERK4, :ESERK5, :RKC]
  @eval begin
    struct $Alg{E} <: OrdinaryDiffEqAdaptiveAlgorithm
      eigen_est::E
    end
    $Alg(;eigen_est=nothing) = $Alg(eigen_est)
  end
end
struct SERK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
  controller::Symbol
  eigen_est::E
end
SERK2(;controller=:PI,eigen_est=nothing) = SERK2(controller,eigen_est)

struct IRKC{CS,AD,F,F2,FDT,K,T,E} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  controller::Symbol
  eigen_est::E
end
IRKC(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),κ=nothing,tol=nothing,
                 extrapolant=:linear,controller = :Standard,eigen_est=nothing) =
  IRKC{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),typeof(κ),typeof(tol),typeof(eigen_est)}(
                 linsolve,nlsolve,diff_type,κ,tol,extrapolant,controller,eigen_est)

################################################################################

# Linear Methods

for Alg in [:MagnusMidpoint,:MagnusLeapfrog,:LieEuler]
  @eval struct $Alg <: OrdinaryDiffEqExponentialAlgorithm
    krylov::Bool
    m::Int
    iop::Int
  end
  @eval $Alg(;krylov=false, m=30, iop=0) = $Alg(krylov, m, iop)
end

struct LinearExponential <: OrdinaryDiffEqExponentialAlgorithm
  krylov::Symbol
  m::Int
  iop::Int
end
LinearExponential(;krylov=:off, m=10, iop=0) = LinearExponential(krylov, m, iop)

struct CayleyEuler <: OrdinaryDiffEqAlgorithm end  

################################################################################

# FIRK Methods

struct RadauIIA5{CS,AD,F,FDT,Tol,C1,C2} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  κ::Tol
  maxiters::Int
  fast_convergence_cutoff::C1
  new_W_γdt_cutoff::C2
  controller::Symbol
end
RadauIIA5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,
                          extrapolant=:dense,fast_convergence_cutoff=1//5,new_W_γdt_cutoff=1//5,
                          controller=:Predictive,κ=nothing,maxiters=10,smooth_est=true) =
                          RadauIIA5{chunk_size,autodiff,typeof(linsolve),
                          typeof(diff_type),
                          typeof(κ),typeof(fast_convergence_cutoff),typeof(new_W_γdt_cutoff)}(
                            linsolve,diff_type,smooth_est,extrapolant,κ,maxiters,fast_convergence_cutoff,new_W_γdt_cutoff,controller)

################################################################################

# SDIRK Methods

struct ImplicitEuler{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
ImplicitEuler(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          controller=:PI) =
                          ImplicitEuler{chunk_size,autodiff,typeof(linsolve),
                          typeof(nlsolve),typeof(diff_type)}(linsolve,
                          nlsolve,diff_type,extrapolant,controller)

struct ImplicitMidpoint{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
ImplicitMidpoint(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      ImplicitMidpoint{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct Trapezoid{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
Trapezoid(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                            extrapolant=:linear,
                      controller = :PI) =
                      Trapezoid{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant,controller)

struct TRBDF2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
TRBDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                 smooth_est=true,extrapolant=:linear,
                 controller = :PI) =
TRBDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
      linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct SDIRK2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
SDIRK2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 SDIRK2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct SDIRK22{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
SDIRK22(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                            extrapolant=:linear,
                      controller = :PI) =
                      Trapezoid{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant,controller)


struct SSPSDIRK2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD} # Not adaptive
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
SSPSDIRK2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:constant,
                   controller = :PI) =
 SSPSDIRK2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct Kvaerno3{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Kvaerno3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Kvaerno3{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct KenCarp3{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp3{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct CFNLIRK3{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CFNLIRK3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CFNLIRK3{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct CFNLIRK4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CFNLIRK4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CFNLIRK4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct Cash4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  embedding::Int
  controller::Symbol
end
Cash4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI,embedding=3) =
 Cash4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,embedding,controller)


struct SFSDIRK4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
SFSDIRK4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

                      struct SFSDIRK5{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
                        linsolve::F
                        nlsolve::F2
                        diff_type::FDT
                        extrapolant::Symbol
                      end
SFSDIRK5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK5{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct SFSDIRK6{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end

SFSDIRK6(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK6{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct SFSDIRK7{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end

SFSDIRK7(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK7{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct SFSDIRK8{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
linsolve::F
nlsolve::F2
diff_type::FDT
extrapolant::Symbol
end

SFSDIRK8(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                    linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                    extrapolant=:linear) =
                    SFSDIRK8{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                    linsolve,nlsolve,diff_type,extrapolant)

struct Hairer4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Hairer4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Hairer4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct Hairer42{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Hairer42(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Hairer42{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct Kvaerno4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Kvaerno4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Kvaerno4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct Kvaerno5{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Kvaerno5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Kvaerno5{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct KenCarp4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct KenCarp5{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp5{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

# `smooth_est` is not necessary, as the embedded method is also L-stable
struct ESDIRK54I8L2SA{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
ESDIRK54I8L2SA(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   extrapolant=:linear,controller = :PI) =
 ESDIRK54I8L2SA{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(linsolve,nlsolve,diff_type,extrapolant,controller)

################################################################################

# Rosenbrock Methods

for Alg in [:Rosenbrock23, :Rosenbrock32, :ROS3P, :Rodas3, :ROS34PW1a, :ROS34PW1b, :ROS34PW2, :ROS34PW3, :RosShamp4, :Veldd4, :Velds4, :GRK4T, :GRK4A, :Ros4LStab, :Rodas4, :Rodas42, :Rodas4P, :Rodas5]
  @eval begin
    struct $Alg{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
      linsolve::F
      diff_type::FDT
    end
    $Alg(;chunk_size=0,autodiff=true,diff_type=Val{:forward},linsolve=DEFAULT_LINSOLVE) = $Alg{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)
  end
end

struct GeneralRosenbrock{CS,AD,F,TabType} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  tableau::TabType
  factorization::F
end

GeneralRosenbrock(;chunk_size=0,autodiff=true,
                    factorization=lu!,tableau=ROSENBROCK_DEFAULT_TABLEAU) =
                    GeneralRosenbrock{chunk_size,autodiff,typeof(factorization),typeof(tableau)}(tableau,factorization)

struct RosenbrockW6S4OS{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
RosenbrockW6S4OS(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = RosenbrockW6S4OS{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)
######################################

for Alg in [:LawsonEuler, :NorsettEuler, :ETDRK2, :ETDRK3, :ETDRK4, :HochOst4]
  @eval struct $Alg{FDT} <: OrdinaryDiffEqExponentialAlgorithm
    krylov::Bool
    m::Int
    iop::Int
    autodiff::Bool
    chunksize::Int
    diff_type::FDT
  end
  @eval $Alg(;krylov=false, m=30, iop=0, autodiff=true, chunksize=0,
            diff_type = Val{:forward}) = $Alg(krylov, m, iop, autodiff,
            chunksize, diff_type)
end
const ETD1 = NorsettEuler # alias
for Alg in [:Exprb32, :Exprb43]
  @eval struct $Alg{FDT} <: OrdinaryDiffEqAdaptiveExponentialAlgorithm
    m::Int
    iop::Int
    autodiff::Bool
    chunksize::Int
    diff_type::FDT
  end
  @eval $Alg(;m=30, iop=0, autodiff=true, chunksize=0,
            diff_type = Val{:forward}) = $Alg(m, iop, autodiff, chunksize, diff_type)
end
for Alg in [:Exp4, :EPIRK4s3A, :EPIRK4s3B, :EPIRK5s3, :EXPRB53s3, :EPIRK5P1, :EPIRK5P2]
  @eval struct $Alg{FDT} <: OrdinaryDiffEqExponentialAlgorithm
    adaptive_krylov::Bool
    m::Int
    iop::Int
    autodiff::Bool
    chunksize::Int
    diff_type::FDT
  end
  @eval $Alg(;adaptive_krylov=true, m=30, iop=0, autodiff=true,
              chunksize=0, diff_type = Val{:forward}) =
              $Alg(adaptive_krylov, m, iop, autodiff, chunksize, diff_type)
end
struct SplitEuler <: OrdinaryDiffEqExponentialAlgorithm end
struct ETD2 <: OrdinaryDiffEqExponentialAlgorithm end

#########################################

struct ABDF2{CS,AD,F,F2,FDT,K,T} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
ABDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
      κ=nothing,tol=nothing,linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
      smooth_est=true,extrapolant=:linear,
      controller=:Standard) =
ABDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
      typeof(κ),typeof(tol)}(
      linsolve,nlsolve,diff_type,κ,tol,smooth_est,extrapolant,controller)

#########################################

struct CompositeAlgorithm{T,F} <: OrdinaryDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

################################################################################

struct MEBDF2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
MEBDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:constant) =
                      MEBDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

#################################################

struct PDIRK44{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  threading::Bool
end
PDIRK44(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:constant,threading=true) =
                      PDIRK44{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant,threading)
### Algorithm Groups

const MultistepAlgorithms = Union{IRKN3,IRKN4,
                                  ABDF2,
                                  AB3,AB4,AB5,ABM32,ABM43,ABM54}

const SplitAlgorithms = Union{CNAB2,CNLF2,IRKC,SBDF,
                              KenCarp3,KenCarp4,KenCarp5,CFNLIRK3,CFNLIRK4}


# DAE Specific Algorithms
abstract type DAEAlgorithm{CS,AD} <: DiffEqBase.AbstractDAEAlgorithm end

#=
struct DBDF{CS,AD,F,F2,FDT} <: DAEAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end

DBDF(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
     linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),extrapolant=:linear) =
     DBDF{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
     linsolve,nlsolve,diff_type,extrapolant)
=#

struct DImplicitEuler{CS,AD,F,F2,FDT} <: DAEAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
DImplicitEuler(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          controller=:Standard) =
                          DImplicitEuler{chunk_size,autodiff,typeof(linsolve),
                          typeof(nlsolve),typeof(diff_type)}(linsolve,
                          nlsolve,diff_type,extrapolant,controller)


struct DABDF2{CS,AD,F,F2,FDT} <: DAEAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
DABDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          controller=:Standard) =
                          DABDF2{chunk_size,autodiff,typeof(linsolve),
                          typeof(nlsolve),typeof(diff_type)}(linsolve,
                          nlsolve,diff_type,extrapolant,controller)
