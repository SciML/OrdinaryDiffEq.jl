abstract AbstractODEAlgorithm
abstract OrdinaryDiffEqAlgorithm <: AbstractODEAlgorithm
abstract OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm

@with_kw immutable Euler <: OrdinaryDiffEqAlgorithm
  order::Int=1
end

@with_kw immutable Midpoint <: OrdinaryDiffEqAlgorithm
  order::Int=2
end

@with_kw immutable RK4 <: OrdinaryDiffEqAlgorithm
  order::Int=4
end

@with_kw immutable ExplicitRK <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=4
  adaptiveorder::Int=3
  tableau::ExplicitRKTableau=constructDormandPrince()
end

@with_kw immutable BS3 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=3
  adaptiveorder::Int=2
end

@with_kw immutable BS5 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=5
  adaptiveorder::Int=4
end

@with_kw immutable DP5 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=5
  adaptiveorder::Int=4
end

@with_kw immutable DP5Threaded <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=5
  adaptiveorder::Int=4
end

@with_kw immutable Tsit5 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=5
  adaptiveorder::Int=4
end

@with_kw immutable DP8 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=8
  adaptiveorder::Int=8
end

@with_kw immutable Vern6 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=6
  adaptiveorder::Int=5
end

@with_kw immutable Vern7 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=7
  adaptiveorder::Int=6
end

@with_kw immutable Vern8 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=8
  adaptiveorder::Int=7
end

@with_kw immutable TanYam7 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=7
  adaptiveorder::Int=6
end

@with_kw immutable TsitPap8 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=8
  adaptiveorder::Int=7
end

@with_kw immutable Vern9 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=9
  adaptiveorder::Int=8
end

@with_kw immutable ImplicitEuler <: OrdinaryDiffEqAlgorithm
  order::Int=1
end

@with_kw immutable Trapezoid <: OrdinaryDiffEqAlgorithm
  order::Int=2
end

@with_kw immutable Rosenbrock23 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=2
  adaptiveorder::Int=2
end

@with_kw immutable Rosenbrock32 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=3
  adaptiveorder::Int=2
end

@with_kw immutable Feagin10 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=10
  adaptiveorder::Int=8
end

@with_kw immutable Feagin12 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=12
  adaptiveorder::Int=10
end

@with_kw immutable Feagin14 <: OrdinaryDiffEqAdaptiveAlgorithm
  order::Int=14
  adaptiveorder::Int=12
end

DefaultODEAlgorithm(alg_hints::Symbol...) = DP5()

function isfsal(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded || typeof(alg) <: DP8 || typeof(alg) <: BS3 || typeof(alg) <: BS5 || typeof(alg) <: Tsit5 || typeof(alg) <: Vern6
    return true
  else
    return false
  end
end

function isspecialdense(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded || typeof(alg) <: DP8 || typeof(alg) <: BS5 || typeof(alg) <: Tsit5 || typeof(alg) <: Vern6 || typeof(alg) <: Vern7 || typeof(alg) <: Vern8 || typeof(alg) <: Vern9
    return true
  else
    return false
  end
end

function isimplicit(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: ImplicitEuler || typeof(alg) <: Trapezoid
    return true
  else
    return false
  end
end

abstract ODEInterfaceAlgorithm <: AbstractODEAlgorithm
immutable dopri5 <: ODEInterfaceAlgorithm end
immutable dop853 <: ODEInterfaceAlgorithm end
immutable odex <: ODEInterfaceAlgorithm end
immutable seulex <: ODEInterfaceAlgorithm end
immutable radau <: ODEInterfaceAlgorithm end
immutable radau5 <: ODEInterfaceAlgorithm end

abstract ODEIterAlgorithm <: AbstractODEAlgorithm
immutable feuler <: ODEIterAlgorithm end
immutable rk23 <: ODEIterAlgorithm end
immutable rk45 <: ODEIterAlgorithm end
immutable feh78 <: ODEIterAlgorithm end
immutable ModifiedRosenbrockIntegrator <: ODEIterAlgorithm end
immutable midpoint <: ODEIterAlgorithm end
immutable heun <: ODEIterAlgorithm end
immutable rk4 <: ODEIterAlgorithm end
immutable feh45 <: ODEIterAlgorithm end

abstract ODEJLAlgorithm <: AbstractODEAlgorithm
immutable ode1 <: ODEIterAlgorithm end
immutable ode23 <: ODEIterAlgorithm end
immutable ode45 <: ODEIterAlgorithm end
immutable ode78 <: ODEIterAlgorithm end
immutable ode23s <: ODEIterAlgorithm end
immutable ode2_midpoint <: ODEIterAlgorithm end
immutable ode2_heun <: ODEIterAlgorithm end
immutable ode4 <: ODEIterAlgorithm end
immutable ode45_fe <: ODEIterAlgorithm end

abstract SundialsAlgorithm <: AbstractODEAlgorithm
immutable cvode_BDF <: SundialsAlgorithm end
immutable cvode_Adams <: SundialsAlgorithm end
