isautodifferentiable(alg::OrdinaryDiffEqAlgorithm) = true

isfsal(alg::OrdinaryDiffEqAlgorithm) = false
isfsal(alg::DP5) = true
isfsal(alg::DP5Threaded) = true
isfsal(alg::DP8) = true
isfsal(alg::BS3) = true
isfsal(alg::BS5) = true
isfsal(alg::Tsit5) = true
isfsal(alg::Vern6) = true
isfsal(alg::Rosenbrock23) = true
isfsal(alg::Rosenbrock32) = true
isfsal(alg::Euler) = true
isfsal(alg::SplitEuler) = true
isfsal(alg::SymplecticEuler) = true
isfsal(alg::Midpoint) = true
isfsal(alg::SSPRK22) = true
isfsal(alg::SSPRK33) = true
isfsal(alg::SSPRK104) = true
isfsal(alg::RK4) = true
isfsal(alg::Feagin10) = true
isfsal(alg::Feagin12) = true
isfsal(alg::Feagin14) = true
isfsal(alg::TanYam7) = true
isfsal(alg::TsitPap8) = true
isfsal(alg::Trapezoid) = true
isfsal(alg::ImplicitEuler) = true
isfsal(alg::ExplicitRK) = true
isfsal{MType,VType,fsal}(tab::ExplicitRKTableau{MType,VType,fsal}) = fsal
#isfsal(tab::ImplicitRKTableau) = false
isfsal(alg::CompositeAlgorithm) = true # Every algorithm is assumed FSAL. Good assumption?

isimplicit(alg::OrdinaryDiffEqAlgorithm) = false
isimplicit(alg::ImplicitEuler) = true
isimplicit(alg::Trapezoid) = true

isdtchangeable(alg::OrdinaryDiffEqAlgorithm) = true

ismultistep(alg::OrdinaryDiffEqAlgorithm) = false

isadaptive(alg::OrdinaryDiffEqAlgorithm) = false
isadaptive(alg::OrdinaryDiffEqAdaptiveAlgorithm) = true
isadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = isadaptive(alg.algs[1])

qmin_default(alg::OrdinaryDiffEqAlgorithm) = 1//5
qmin_default(alg::DP8) = 1//3

qmax_default(alg::OrdinaryDiffEqAlgorithm) = 10
qmax_default(alg::DP8) = 6

get_chunksize(alg::OrdinaryDiffEqAlgorithm) = error("This algorithm does not have a chunk size defined.")
get_chunksize{CS,AD}(alg::ImplicitEuler{CS,AD}) = CS
get_chunksize{CS,AD}(alg::Trapezoid{CS,AD}) = CS
get_chunksize{CS,AD}(alg::Rosenbrock23{CS,AD}) = CS
get_chunksize{CS,AD}(alg::Rosenbrock32{CS,AD}) = CS

alg_extrapolates(alg::OrdinaryDiffEqAlgorithm) = false
alg_extrapolates(alg::ImplicitEuler) = true
alg_extrapolates(alg::Trapezoid) = true

alg_autodiff(alg::OrdinaryDiffEqAlgorithm) = error("This algorithm does not have an autodifferentiation option defined.")
alg_autodiff{CS,AD}(alg::ImplicitEuler{CS,AD}) = AD
alg_autodiff{CS,AD}(alg::Trapezoid{CS,AD}) = AD
alg_autodiff{CS,AD}(alg::Rosenbrock23{CS,AD}) = AD
alg_autodiff{CS,AD}(alg::Rosenbrock32{CS,AD}) = AD

alg_order(alg::OrdinaryDiffEqAlgorithm) = error("Order is not defined for this algorithm")
alg_adaptive_order(alg::OrdinaryDiffEqAdaptiveAlgorithm) = error("Algorithm is adaptive with no order")

alg_order(alg::Discrete) = 0
alg_order(alg::Euler) = 1
alg_order(alg::SplitEuler) = 1
alg_order(alg::SymplecticEuler) = 1
alg_order(alg::Midpoint) = 2
alg_order(alg::SSPRK22) = 2
alg_order(alg::SSPRK33) = 3
alg_order(alg::SSPRK104) = 4
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

alg_order(alg::CompositeAlgorithm) = alg_order(alg.algs[1])

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

beta2_default(alg::OrdinaryDiffEqAlgorithm) = 2//(5alg_order(alg))
beta2_default(alg::Discrete) = 0
beta2_default(alg::DP8) = 0//1
beta2_default(alg::DP5) = 4//100
beta2_default(alg::DP5Threaded) = 4//100

beta1_default(alg::OrdinaryDiffEqAlgorithm,beta2) = 7//(10alg_order(alg))
beta1_default(alg::Discrete,beta2) = 0
beta1_default(alg::DP8,beta2) = typeof(beta2)(1//alg_order(alg)) - beta2/5
beta1_default(alg::DP5,beta2) = typeof(beta2)(1//alg_order(alg)) - 3beta2/4


discrete_apply_map{apply_map,scale_by_time}(alg::Discrete{apply_map,scale_by_time}) = apply_map
discrete_scale_by_time{apply_map,scale_by_time}(alg::Discrete{apply_map,scale_by_time}) = scale_by_time
