abstract type OrdinaryDiffEqAlgorithm <: AbstractODEAlgorithm end
abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

struct Discrete{apply_map,scale_by_time} <: OrdinaryDiffEqAlgorithm end

Base.@pure Discrete(;apply_map=false,scale_by_time=false) = Discrete{apply_map,scale_by_time}()
Base.@pure FunctionMap(;scale_by_time=false) = Discrete{true,scale_by_time}()

###############################################################################

# RK methods

@with_kw struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
  tableau::TabType=ODE_DEFAULT_TABLEAU
end

struct Euler <: OrdinaryDiffEqAlgorithm end
struct Midpoint <: OrdinaryDiffEqAlgorithm end
struct RK4 <: OrdinaryDiffEqAlgorithm end
struct SSPRK22 <: OrdinaryDiffEqAlgorithm end
struct SSPRK33 <: OrdinaryDiffEqAlgorithm end
struct SSPRK104 <: OrdinaryDiffEqAlgorithm end
struct SSPRK432 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct BS3 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct BS5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP5Threaded <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Tsit5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct DP8 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Vern6 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Vern7 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Vern8 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct TanYam7 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct TsitPap8 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Vern9 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Feagin10 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Feagin12 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Feagin14 <: OrdinaryDiffEqAdaptiveAlgorithm end

################################################################################

# Symplectic methods

#struct Verlet <: OrdinaryDiffEqAlgorithm end
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

struct Nystrom4 <: OrdinaryDiffEqAlgorithm end
struct Nystrom4VelocityIndependent <: OrdinaryDiffEqAlgorithm end
struct Nystrom5VelocityIndependent <: OrdinaryDiffEqAlgorithm end

################################################################################

# Generic implicit methods

struct GenericImplicitEuler{F} <: OrdinaryDiffEqAdaptiveAlgorithm
  nlsolve::F
  extrapolant::Symbol
end
Base.@pure GenericImplicitEuler(;
            nlsolve=NLSOLVEJL_SETUP(),extrapolant=:constant) =
            GenericImplicitEuler{typeof(nlsolve)}(nlsolve,extrapolant)

struct GenericTrapezoid{F} <: OrdinaryDiffEqAdaptiveAlgorithm
  nlsolve::F
  extrapolant::Symbol
end
Base.@pure GenericTrapezoid(;
            nlsolve=NLSOLVEJL_SETUP(),extrapolant=:constant) =
            GenericTrapezoid{typeof(nlsolve)}(nlsolve,extrapolant)

################################################################################

# Linear Methods

struct LinearImplicitEuler{F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
end
Base.@pure LinearImplicitEuler(;linsolve=DEFAULT_LINSOLVE) = LinearImplicitEuler{typeof(linsolve)}(linsolve)

################################################################################

# Implicit RK Methods

struct ImplicitEuler{CS,AD,F,K,T} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
  κ::K
  tol::T
  extrapolant::Symbol
end
Base.@pure ImplicitEuler(;chunk_size=0,autodiff=true,diff_type=:central,
                          linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                          extrapolant=:constant) = ImplicitEuler{chunk_size,autodiff,typeof(linsolve),
                          typeof(κ),typeof(tol)}(
                          linsolve,diff_type,κ,tol,extrapolant)

struct Trapezoid{CS,AD,F,K,T} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
  κ::K
  tol::T
  extrapolant::Symbol
end
Base.@pure Trapezoid(;chunk_size=0,autodiff=true,diff_type=:central,
                      linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                      extrapolant=:constant) =
                      Trapezoid{chunk_size,autodiff,typeof(linsolve),
                      typeof(κ),typeof(tol)}(
                      linsolve,diff_type,κ,tol,extrapolant)

struct TRBDF2{CS,AD,F,K,T} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
end
Base.@pure TRBDF2(;chunk_size=0,autodiff=true,diff_type=:central,
                   linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                   smooth_est=true,extrapolant=:constant) =
 TRBDF2{chunk_size,autodiff,typeof(linsolve),typeof(κ),typeof(tol)}(
        linsolve,diff_type,κ,tol,smooth_est,extrapolant)

################################################################################

# Rosenbrock Methods

struct Rosenbrock23{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rosenbrock23(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rosenbrock23{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Rosenbrock32{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rosenbrock32(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rosenbrock32{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct ROS3P{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure ROS3P(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = ROS3P{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Rodas3{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rodas3(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rodas3{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct RosShamp4{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure RosShamp4(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = RosShamp4{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Veldd4{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Veldd4(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Veldd4{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Velds4{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Velds4(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Velds4{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct GRK4T{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure GRK4T(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = GRK4T{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct GRK4A{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure GRK4A(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = GRK4A{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Ros4LStab{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Ros4LStab(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Ros4LStab{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Rodas4{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rodas4(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rodas4{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Rodas42{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rodas42(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rodas42{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Rodas4P{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rodas4P(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rodas4P{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct Rodas5{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  linsolve::F
  diff_type::Symbol
end
Base.@pure Rodas5(;chunk_size=0,autodiff=true,diff_type=:central,linsolve=DEFAULT_LINSOLVE) = Rodas5{chunk_size,autodiff,typeof(linsolve)}(linsolve,diff_type)

struct GeneralRosenbrock{CS,AD,F,TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
  tableau::TabType
  factorization::F
end

Base.@pure GeneralRosenbrock(;chunk_size=0,autodiff=true,
                    factorization=lufact!,tableau=ROSENBROCK_DEFAULT_TABLEAU) =
                    GeneralRosenbrock{chunk_size,autodiff,typeof(factorization),typeof(tableau)}(tableau,factorization)

######################################

struct IIF1{F} <: OrdinaryDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF1(;nlsolve=NLSOLVEJL_SETUP()) = IIF1{typeof(nlsolve)}(nlsolve)

struct IIF2{F} <: OrdinaryDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF2(;nlsolve=NLSOLVEJL_SETUP()) = IIF2{typeof(nlsolve)}(nlsolve)

struct LawsonEuler <: OrdinaryDiffEqAlgorithm end
struct NorsettEuler <: OrdinaryDiffEqAlgorithm end
struct SplitEuler <: OrdinaryDiffEqAlgorithm end

struct CompositeAlgorithm{T,F} <: OrdinaryDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end
