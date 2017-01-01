abstract OrdinaryDiffEqAlgorithm <: AbstractODEAlgorithm
abstract OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm

immutable Euler <: OrdinaryDiffEqAlgorithm end
immutable Midpoint <: OrdinaryDiffEqAlgorithm end
immutable RK4 <: OrdinaryDiffEqAlgorithm end

@with_kw immutable ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
  tableau::TabType=ODE_DEFAULT_TABLEAU
end

immutable BS3 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable BS5 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable DP5 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable DP5Threaded <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Tsit5 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable DP8 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Vern6 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Vern7 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Vern8 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable TanYam7 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable TsitPap8 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Vern9 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Feagin10 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Feagin12 <: OrdinaryDiffEqAdaptiveAlgorithm end
immutable Feagin14 <: OrdinaryDiffEqAdaptiveAlgorithm end

immutable ImplicitEuler{CS,AD} <: OrdinaryDiffEqAlgorithm end
Base.@pure ImplicitEuler(;chunk_size=1,autodiff=true) = ImplicitEuler{chunk_size,autodiff}()

immutable Trapezoid{CS,AD} <: OrdinaryDiffEqAlgorithm end
Base.@pure Trapezoid(;chunk_size=1,autodiff=true) = Trapezoid{chunk_size,autodiff}()

immutable Rosenbrock23{CS,AD} <: OrdinaryDiffEqAdaptiveAlgorithm end
Base.@pure Rosenbrock23(;chunk_size=1,autodiff=true) = Rosenbrock23{chunk_size,autodiff}()

immutable Rosenbrock32{CS,AD} <: OrdinaryDiffEqAdaptiveAlgorithm end
Base.@pure Rosenbrock32(;chunk_size=1,autodiff=true) = Rosenbrock32{chunk_size,autodiff}()
