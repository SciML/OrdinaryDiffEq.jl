function isfsal(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded || typeof(alg) <: DP8 || typeof(alg) <: BS3 || typeof(alg) <: BS5 || typeof(alg) <: Tsit5 || typeof(alg) <: Vern6 || typeof(alg) <: Rosenbrock23 || typeof(alg) <: Rosenbrock32 || typeof(alg)<:Euler || typeof(alg) <: Midpoint || typeof(alg) <: RK4 || typeof(alg) <: Feagin10 || typeof(alg) <: Feagin12 || typeof(alg) <: Feagin14
    return true
  else
    return false
  end
end

isfsal(alg::ExplicitRK) = isfsal(alg.tableau)
isfsal{MType,VType,fsal}(tab::ExplicitRKTableau{MType,VType,fsal}) = fsal
isfsal(tab::ImplicitRKTableau) = false

function isspecialdense(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded || typeof(alg) <: DP8 || typeof(alg) <: BS5 || typeof(alg) <: Tsit5 || typeof(alg) <: Vern6 || typeof(alg) <: Vern7 || typeof(alg) <: Vern8 || typeof(alg) <: Vern9 ||
  typeof(alg) <: Rosenbrock23 || typeof(alg) <: Rosenbrock32
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

function get_kseltype(alg::OrdinaryDiffEqAlgorithm,prob)
  rateType = typeof(prob.u0/zero(prob.tspan[1]))
  isspecialdense(alg) ? ksEltype = Vector{rateType} : ksEltype = rateType
end

qmin_default(alg::OrdinaryDiffEqAlgorithm) = 0.2
qmin_default(alg::DP8) = 0.333

qmax_default(alg::OrdinaryDiffEqAlgorithm) = 10.0
qmax_default(alg::DP8) = 6.0

get_chunksize(alg::OrdinaryDiffEqAlgorithm) = error("This algorithm does not have a chunk size defined.")
get_chunksize{CS,AD}(alg::ImplicitEuler{CS,AD}) = CS
get_chunksize{CS,AD}(alg::Trapezoid{CS,AD}) = CS
get_chunksize{CS,AD}(alg::Rosenbrock23{CS,AD}) = CS
get_chunksize{CS,AD}(alg::Rosenbrock32{CS,AD}) = CS

alg_autodiff(alg::OrdinaryDiffEqAlgorithm) = error("This algorithm does not have an autodifferentiation option defined.")
alg_autodiff{CS,AD}(alg::ImplicitEuler{CS,AD}) = AD
alg_autodiff{CS,AD}(alg::Trapezoid{CS,AD}) = AD
alg_autodiff{CS,AD}(alg::Rosenbrock23{CS,AD}) = AD
alg_autodiff{CS,AD}(alg::Rosenbrock32{CS,AD}) = AD


alg_order(alg::OrdinaryDiffEqAlgorithm) = error("Order is not defined for this algorithm")
alg_adaptive_order(alg::OrdinaryDiffEqAdaptiveAlgorithm) = error("Algorithm is adaptive with no order")

alg_order(alg::Euler) = 1
alg_order(alg::Midpoint) = 2
alg_order(alg::RK4) = 4
alg_order(alg::ExplicitRK) = alg.tableau.order
alg_order(alg::BS3) = 3
alg_order(alg::BS5) = 5
alg_order(alg::DP5) = 5
alg_order(alg::DP5Threaded) = 5
alg_order(alg::Tsit5) = 5
alg_order(alg::DP8) = 8
alg_order(alg::Vern6) = 6
alg_order(alg::Vern7) = 7
alg_order(alg::Vern8) = 8
alg_order(alg::Vern9) = 9
alg_order(alg::TanYam7) = 7
alg_order(alg::TsitPap8) = 8
alg_order(alg::ImplicitEuler) = 1
alg_order(alg::Trapezoid) = 2
alg_order(alg::Rosenbrock23) = 2
alg_order(alg::Rosenbrock32) = 3
alg_order(alg::Feagin10) = 10
alg_order(alg::Feagin12) = 12
alg_order(alg::Feagin14) = 14

alg_adaptive_order(alg::ExplicitRK) = alg.tableau.adaptiveorder
alg_adaptive_order(alg::BS3) = 2
alg_adaptive_order(alg::BS5) = 4
alg_adaptive_order(alg::DP5) = 4
alg_adaptive_order(alg::DP5Threaded) = 4
alg_adaptive_order(alg::Tsit5) = 4
alg_adaptive_order(alg::DP8) = 6
alg_adaptive_order(alg::Vern6) = 5
alg_adaptive_order(alg::Vern7) = 6
alg_adaptive_order(alg::Vern8) = 7
alg_adaptive_order(alg::Vern9) = 8
alg_adaptive_order(alg::TanYam7) = 6
alg_adaptive_order(alg::TsitPap8) = 7
alg_adaptive_order(alg::Rosenbrock23) = 3
alg_adaptive_order(alg::Rosenbrock32) = 2
alg_adaptive_order(alg::Feagin10) = 8
alg_adaptive_order(alg::Feagin12) = 10
alg_adaptive_order(alg::Feagin14) = 12

beta2_default(alg::OrdinaryDiffEqAlgorithm) = 0.4/alg_order(alg)
beta2_default(alg::DP8) = 0.0
beta2_default(alg::DP5) = 0.04
beta2_default(alg::DP5Threaded) = 0.04

beta1_default(alg::OrdinaryDiffEqAlgorithm,beta2) = 0.7/alg_order(alg)
beta1_default(alg::DP8,beta2) = 1/alg_order(alg) - 0.2beta2
beta1_default(alg::DP5,beta2) = 1/alg_order(alg) - 0.75beta2
