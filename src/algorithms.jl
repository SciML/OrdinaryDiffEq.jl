@compat abstract type OrdinaryDiffEqAlgorithm <: AbstractODEAlgorithm end
@compat abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
@compat abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

immutable Discrete{apply_map,scale_by_time} <: OrdinaryDiffEqAlgorithm end

Base.@pure Discrete(;apply_map=false,scale_by_time=false) = Discrete{apply_map,scale_by_time}()
Base.@pure FunctionMap(;scale_by_time=false) = Discrete{true,scale_by_time}()
immutable Euler <: OrdinaryDiffEqAlgorithm end
immutable Midpoint <: OrdinaryDiffEqAlgorithm end
immutable RK4 <: OrdinaryDiffEqAlgorithm end
immutable SSPRK22 <: OrdinaryDiffEqAlgorithm end

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

immutable ImplicitEuler{CS,AD,F} <: OrdinaryDiffEqAlgorithm
  factorization::F
end
Base.@pure ImplicitEuler(;chunk_size=0,autodiff=true,factorization=lufact!) = ImplicitEuler{chunk_size,autodiff,typeof(factorization)}(factorization)

immutable Trapezoid{CS,AD,F} <: OrdinaryDiffEqAlgorithm
  factorization::F
end
Base.@pure Trapezoid(;chunk_size=0,autodiff=true,factorization=lufact!) = Trapezoid{chunk_size,autodiff,typeof(factorization)}(factorization)

immutable Rosenbrock23{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  factorization::F
end
Base.@pure Rosenbrock23(;chunk_size=0,autodiff=true,factorization=lufact!) = Rosenbrock23{chunk_size,autodiff,typeof(factorization)}(factorization)

immutable Rosenbrock32{CS,AD,F} <: OrdinaryDiffEqAdaptiveAlgorithm
  factorization::F
end
Base.@pure Rosenbrock32(;chunk_size=0,autodiff=true,factorization=lufact!) = Rosenbrock32{chunk_size,autodiff,typeof(factorization)}(factorization)

immutable CompositeAlgorithm{T,F} <: OrdinaryDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end
