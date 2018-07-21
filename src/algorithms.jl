abstract type OrdinaryDiffEqAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

abstract type OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end
abstract type OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end

abstract type OrdinaryDiffEqImplicitAlgorithm{CS,AD} <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqNewtonAlgorithm{CS,AD} <:  OrdinaryDiffEqImplicitAlgorithm{CS,AD} end
abstract type OrdinaryDiffEqRosenbrockAlgorithm{CS,AD} <:  OrdinaryDiffEqImplicitAlgorithm{CS,AD} end

abstract type OrdinaryDiffEqExponentialAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqAdaptiveExponentialAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end

abstract type OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end

struct FunctionMap{scale_by_time} <: OrdinaryDiffEqAlgorithm end
Base.@pure function Discrete(;apply_map=false,scale_by_time=false)
    @warn("Discrete is deprecated. Use FunctionMap instead.")
    FunctionMap{scale_by_time}()
end
Base.@pure FunctionMap(;scale_by_time=false) = FunctionMap{scale_by_time}()

###############################################################################

# RK methods

@with_kw struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
  tableau::TabType=ODE_DEFAULT_TABLEAU
end

@inline trivial_limiter!(u, f, t) = nothing

struct Euler <: OrdinaryDiffEqAlgorithm end
struct Heun <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Ralston <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Midpoint <: OrdinaryDiffEqAdaptiveAlgorithm end
struct RK4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Anas5{T} <: OrdinaryDiffEqAlgorithm
  w::T		
end
Base.@pure Anas5(; w=1) = Anas5(w)

struct OwrenZen3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct OwrenZen4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct OwrenZen5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CarpenterKennedy2N54 <: OrdinaryDiffEqAlgorithm end
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
SSPRK53(stage_limiter! = trivial_limiter!) = SSPRK53(stage_limiter!, trivial_limiter!)
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
struct BS3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP5Threaded <: OrdinaryDiffEqAdaptiveAlgorithm end
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

# IMEXEuler

struct IMEXEuler{CS,AD,F,FDT,K,T,T2} <: OrdinaryDiffEqImplicitAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure IMEXEuler(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3) =
                      IMEXEuler{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound)}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound)


# IMEX Multistep methods

struct CNAB2{CS,AD,F,FDT,K,T,T2} <: OrdinaryDiffEqImplicitAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure CNAB2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3) =
                      CNAB2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound)}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound)

struct CNLF2{CS,AD,F,FDT,K,T,T2} <: OrdinaryDiffEqImplicitAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure CNLF2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3) =
                      CNLF2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound)}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound)

struct QNDF1{CS,AD,F,FDT,K,T,T2,κType,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
  kappa::κType
end
Base.@pure QNDF1(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3,kappa = -0.1850,
                      controller = :Predictive) =
                      QNDF1{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound),typeof(kappa),controller}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound,kappa)

Base.@pure QBDF1(;kwargs...) = QNDF1(;kappa=0,kwargs...)

struct QNDF2{CS,AD,F,FDT,K,T,T2,κType,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
  kappa::κType
end
Base.@pure QNDF2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3,kappa = -1//9,
                      controller = :Predictive) =
                      QNDF2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound),typeof(kappa),controller}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound,kappa)

Base.@pure QBDF2(;kwargs...) = QNDF2(;kappa=0,kwargs...)


struct QNDF{CS,AD,F,FDT,K,T,T2,κType,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
  kappa::κType
end
Base.@pure QNDF(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3,kappa = tuple(-0.1850,-1//9,-0.0823,-0.0415,0),
                      controller = :Predictive) =
                      QNDF{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound),typeof(kappa),controller}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound,kappa)

Base.@pure QBDF(;kwargs...) = QNDF(;kappa=tuple(0,0,0,0,0),kwargs...)

# Adams/BDF methods in Nordsieck forms
struct AN5   <: OrdinaryDiffEqAdaptiveAlgorithm end
struct JVODE{bType,aType} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
  algorithm::Symbol
  bias1::bType
  bias2::bType
  bias3::bType
  addon::aType
end

Base.@pure JVODE(algorithm=:Adams;bias1=6, bias2=6,bias3=10,
                 addon=1//10^6) = JVODE(algorithm,bias1,bias2,bias3,addon)
Base.@pure JVODE_Adams(;kwargs...) = JVODE(:Adams;kwargs...)

# ROCK methods
struct ROCK2 <: OrdinaryDiffEqAdaptiveAlgorithm end

################################################################################

# Generic implicit methods

struct GenericImplicitEuler{F} <: OrdinaryDiffEqAdaptiveAlgorithm
  nlsolve::F
  extrapolant::Symbol
end
Base.@pure GenericImplicitEuler(;
            nlsolve=NLSOLVEJL_SETUP(),extrapolant=:linear) =
            GenericImplicitEuler{typeof(nlsolve)}(nlsolve,extrapolant)

struct GenericTrapezoid{F} <: OrdinaryDiffEqAdaptiveAlgorithm
  nlsolve::F
  extrapolant::Symbol
end
Base.@pure GenericTrapezoid(;
            nlsolve=NLSOLVEJL_SETUP(),extrapolant=:linear) =
            GenericTrapezoid{typeof(nlsolve)}(nlsolve,extrapolant)

################################################################################

# Linear Methods

struct LinearImplicitEuler{F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
end
Base.@pure LinearImplicitEuler(;linsolve=DEFAULT_LINSOLVE) = LinearImplicitEuler{typeof(linsolve)}(linsolve)

struct MidpointSplitting <: OrdinaryDiffEqAlgorithm end

################################################################################

# SDIRK Methods

struct ImplicitEuler{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure ImplicitEuler(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                          extrapolant=:constant,min_newton_iter=1,
                          max_newton_iter=7,new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEuler{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                          typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
                          linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                          max_newton_iter,new_jac_conv_bound)

struct ImplicitMidpoint{CS,AD,F,FDT,K,T,T2} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure ImplicitMidpoint(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:constant,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3) =
                      ImplicitMidpoint{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound)}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound)

struct Trapezoid{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure Trapezoid(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:linear,min_newton_iter=1,
                      max_newton_iter=7,new_jac_conv_bound = 1e-3,
                      controller = :PI) =
                      Trapezoid{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
                      typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
                      linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                      max_newton_iter,new_jac_conv_bound)

struct TRBDF2{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure TRBDF2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                 linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                 smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                 max_newton_iter=7,new_jac_conv_bound = 1e-3,
                 controller = :Predictive) =
TRBDF2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
      typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
      linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
      max_newton_iter,new_jac_conv_bound)

struct SDIRK2{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure SDIRK2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 SDIRK2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct SSPSDIRK2{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD} # Not adaptive
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure SSPSDIRK2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:constant,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 SSPSDIRK2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct Kvaerno3{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure Kvaerno3(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 Kvaerno3{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct KenCarp3{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure KenCarp3(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 KenCarp3{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)


struct Cash4{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
  embedding::Int
end
Base.@pure Cash4(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive,embedding=3) =
 Cash4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound,embedding)

struct Hairer4{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure Hairer4(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 Hairer4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct Hairer42{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure Hairer42(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 Hairer42{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct Kvaerno4{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure Kvaerno4(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 Kvaerno4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct Kvaerno5{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure Kvaerno5(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 Kvaerno5{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct KenCarp4{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure KenCarp4(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 KenCarp4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

struct KenCarp5{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure KenCarp5(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                   max_newton_iter=7,new_jac_conv_bound = 1e-3,
                   controller = :Predictive) =
 KenCarp5{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
        max_newton_iter,new_jac_conv_bound)

################################################################################

# Rosenbrock Methods

struct Rosenbrock23{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rosenbrock23(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rosenbrock23{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Rosenbrock32{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rosenbrock32(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rosenbrock32{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct ROS3P{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure ROS3P(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = ROS3P{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Rodas3{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rodas3(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rodas3{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct RosShamp4{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure RosShamp4(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = RosShamp4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Veldd4{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Veldd4(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Veldd4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Velds4{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Velds4(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Velds4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct GRK4T{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure GRK4T(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = GRK4T{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct GRK4A{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure GRK4A(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = GRK4A{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Ros4LStab{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Ros4LStab(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Ros4LStab{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Rodas4{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rodas4(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rodas4{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Rodas42{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rodas42(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rodas42{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Rodas4P{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rodas4P(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rodas4P{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct Rodas5{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
Base.@pure Rodas5(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = Rodas5{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)

struct GeneralRosenbrock{CS,AD,F,TabType} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  tableau::TabType
  factorization::F
end

Base.@pure GeneralRosenbrock(;chunk_size=0,autodiff=true,
                    factorization=lu!,tableau=ROSENBROCK_DEFAULT_TABLEAU) =
                    GeneralRosenbrock{chunk_size,autodiff,typeof(factorization),typeof(tableau)}(tableau,factorization)

######################################

struct GenericIIF1{F} <: OrdinaryDiffEqExponentialAlgorithm
  nlsolve::F
end
Base.@pure GenericIIF1(;nlsolve=NLSOLVEJL_SETUP()) = GenericIIF1{typeof(nlsolve)}(nlsolve)

struct GenericIIF2{F} <: OrdinaryDiffEqExponentialAlgorithm
  nlsolve::F
end
Base.@pure GenericIIF2(;nlsolve=NLSOLVEJL_SETUP()) = GenericIIF2{typeof(nlsolve)}(nlsolve)

for Alg in [:LawsonEuler, :NorsettEuler, :ETDRK2, :ETDRK3, :ETDRK4, :HochOst4]
  @eval struct $Alg <: OrdinaryDiffEqExponentialAlgorithm
    krylov::Bool
    m::Int
    iop::Int
  end
  @eval Base.@pure $Alg(;krylov=false, m=30, iop=0) = $Alg(krylov, m, iop)
end
ETD1 = NorsettEuler # alias
for Alg in [:Exp4, :EPIRK4s3A, :EPIRK4s3B, :EPIRK5s3, :EXPRB53s3, :EPIRK5P1, :EPIRK5P2]
  @eval struct $Alg <: OrdinaryDiffEqExponentialAlgorithm
    m::Int
    iop::Int
  end
  @eval Base.@pure $Alg(;m=30, iop=0) = $Alg(m, iop)
end
struct SplitEuler <: OrdinaryDiffEqExponentialAlgorithm end
struct ETD2 <: OrdinaryDiffEqExponentialAlgorithm end

#########################################

struct ABDF2{CS,AD,F,FDT,K,T,T2,Controller} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure ABDF2(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                 linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                 smooth_est=true,extrapolant=:linear,min_newton_iter=1,
                 max_newton_iter=7,new_jac_conv_bound = 1e-3,
                 controller = :Predictive) =
ABDF2{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
      typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
      linsolve,diff_type,κ,tol,smooth_est,extrapolant,min_newton_iter,
      max_newton_iter,new_jac_conv_bound)

#########################################

struct CompositeAlgorithm{T,F} <: OrdinaryDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

################################################################################

### Algorithm Groups

const MassMatrixAlgorithms = Union{OrdinaryDiffEqRosenbrockAlgorithm,
                                   OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                                   ImplicitEuler,ImplicitMidpoint}

const MultistepAlgorithms = Union{IRKN3,IRKN4,
                                  ABDF2,
                                  AB3,AB4,AB5,ABM32,ABM43,ABM54}
