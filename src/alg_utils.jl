isautodifferentiable(alg::OrdinaryDiffEqAlgorithm) = true

isfsal(alg::OrdinaryDiffEqAlgorithm) = true
isfsal(tab::DiffEqBase.ExplicitRKTableau{MType,VType,fsal}) where {MType,VType,fsal} = fsal
# isfsal(alg::CompositeAlgorithm) = isfsal(alg.algs[alg.current])
isfsal(alg::FunctionMap) = false
isfsal(alg::Rodas4) = false
isfsal(alg::Rodas42) = false
isfsal(alg::Rodas4P) = false
isfsal(alg::Vern7) = false
isfsal(alg::Vern8) = false
isfsal(alg::Vern9) = false
get_current_isfsal(alg, cache) = isfsal(alg)
get_current_isfsal(alg::CompositeAlgorithm, cache) = isfsal(alg.algs[cache.current])

fsal_typeof(alg::OrdinaryDiffEqAlgorithm,rate_prototype) = typeof(rate_prototype)
fsal_typeof(alg::ETD2,rate_prototype) = ETD2Fsal{typeof(rate_prototype)}
function fsal_typeof(alg::CompositeAlgorithm,rate_prototype)
  fsal = unique(map(x->fsal_typeof(x,rate_prototype), alg.algs))
  @assert length(fsal) == 1 "`fsal_typeof` must be consistent"
  return fsal[1]
end

isimplicit(alg::OrdinaryDiffEqAlgorithm) = false
isimplicit(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
isimplicit(alg::OrdinaryDiffEqImplicitAlgorithm) = true
isimplicit(alg::CompositeAlgorithm) = any(isimplicit.(alg.algs))

isdtchangeable(alg::OrdinaryDiffEqAlgorithm) = true
isdtchangeable(alg::CompositeAlgorithm) = all(isdtchangeable.(alg.algs))
isdtchangeable(alg::GenericIIF1) = false
isdtchangeable(alg::GenericIIF2) = false
isdtchangeable(alg::LawsonEuler) = false
isdtchangeable(alg::NorsettEuler) = false
isdtchangeable(alg::ETD2) = false

ismultistep(alg::OrdinaryDiffEqAlgorithm) = false
ismultistep(alg::CompositeAlgorithm) = any(ismultistep.(alg.algs))
ismultistep(alg::ETD2) = true

isadaptive(alg::OrdinaryDiffEqAlgorithm) = false
isadaptive(alg::OrdinaryDiffEqAdaptiveAlgorithm) = true
isadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = all(isadaptive.(alg.algs))

qmin_default(alg::OrdinaryDiffEqAlgorithm) = 1//5
qmin_default(alg::CompositeAlgorithm) = maximum(qmin_default.(alg.algs))
qmin_default(alg::DP8) = 1//3

qmax_default(alg::OrdinaryDiffEqAlgorithm) = 10
qmax_default(alg::CompositeAlgorithm) = minimum(qmax_default.(alg.algs))
qmax_default(alg::DP8) = 6

get_chunksize(alg::OrdinaryDiffEqAlgorithm) = error("This algorithm does not have a chunk size defined.")
get_chunksize(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD}) where {CS,AD} = CS
get_chunksize(alg::OrdinaryDiffEqImplicitAlgorithm{CS,AD}) where {CS,AD} = CS
# get_chunksize(alg::CompositeAlgorithm) = get_chunksize(alg.algs[alg.current_alg])

alg_autodiff(alg::OrdinaryDiffEqAlgorithm) = error("This algorithm does not have an autodifferentiation option defined.")
alg_autodiff(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD}) where {CS,AD} = AD
alg_autodiff(alg::OrdinaryDiffEqImplicitAlgorithm{CS,AD}) where {CS,AD} = AD
# alg_autodiff(alg::CompositeAlgorithm) = alg_autodiff(alg.algs[alg.current_alg])
get_current_alg_autodiff(alg, cache) = alg_autodiff(alg)
get_current_alg_autodiff(alg::CompositeAlgorithm, cache) = alg_autodiff(alg.algs[cache.current])

alg_extrapolates(alg::OrdinaryDiffEqAlgorithm) = false
alg_extrapolates(alg::CompositeAlgorithm) = any(alg_extrapolates.(alg.algs))
alg_extrapolates(alg::GenericImplicitEuler) = true
alg_extrapolates(alg::GenericTrapezoid) = true
alg_extrapolates(alg::ImplicitEuler) = true
alg_extrapolates(alg::IMEXEuler) = true
alg_extrapolates(alg::LinearImplicitEuler) = true
alg_extrapolates(alg::Trapezoid) = true
alg_extrapolates(alg::ImplicitMidpoint) = true
alg_extrapolates(alg::TRBDF2) = true
alg_extrapolates(alg::SSPSDIRK2) = true
alg_extrapolates(alg::SDIRK2) = true
alg_extrapolates(alg::Kvaerno3) = true
alg_extrapolates(alg::Kvaerno4) = true
alg_extrapolates(alg::Kvaerno5) = true
alg_extrapolates(alg::KenCarp3) = true
alg_extrapolates(alg::KenCarp4) = true
alg_extrapolates(alg::KenCarp5) = true
alg_extrapolates(alg::Cash4) = true
alg_extrapolates(alg::Hairer4) = true
alg_extrapolates(alg::Hairer42) = true
alg_extrapolates(alg::IRKN4) = true
alg_extrapolates(alg::IRKN3) = true
alg_extrapolates(alg::ABDF2) = true

alg_order(alg::OrdinaryDiffEqAlgorithm) = error("Order is not defined for this algorithm")
get_current_alg_order(alg::OrdinaryDiffEqAlgorithm,cache) = alg_order(alg)
get_current_alg_order(alg::CompositeAlgorithm,cache) = alg_order(alg.algs[cache.current])

get_current_alg_order(alg::OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm,cache) = cache.order
get_current_adaptive_order(alg::OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm,cache) = cache.order
get_current_alg_order(alg::JVODE,cache) = get_current_adaptive_order(alg,cache)

alg_adaptive_order(alg::OrdinaryDiffEqAdaptiveAlgorithm) = error("Algorithm is adaptive with no order")
get_current_adaptive_order(alg::OrdinaryDiffEqAlgorithm,cache) = alg_adaptive_order(alg)
get_current_adaptive_order(alg::CompositeAlgorithm,cache) = alg_adaptive_order(alg.algs[cache.current])

alg_order(alg::FunctionMap) = 0
alg_order(alg::Euler) = 1
alg_order(alg::Heun) = 2
alg_order(alg::Ralston) = 2
alg_order(alg::LawsonEuler) = 1
alg_order(alg::NorsettEuler) = 1
alg_order(alg::ETDRK2) = 2
alg_order(alg::ETDRK3) = 3
alg_order(alg::ETDRK4) = 4
alg_order(alg::HochOst4) = 4
alg_order(alg::Exp4) = 4
alg_order(alg::EPIRK4s3A) = 4
alg_order(alg::EPIRK4s3B) = 4
alg_order(alg::EPIRK5s3) = 5
alg_order(alg::EPIRK5P1) = 5
alg_order(alg::EPIRK5P2) = 5
alg_order(alg::EXPRB53s3) = 5
alg_order(alg::SplitEuler) = 1
alg_order(alg::ETD2) = 2
alg_order(alg::Anas5) = 5

alg_order(alg::SymplecticEuler) = 1
alg_order(alg::VelocityVerlet) = 2
alg_order(alg::VerletLeapfrog) = 2
alg_order(alg::PseudoVerletLeapfrog) = 2
alg_order(alg::McAte2) = 2
alg_order(alg::Ruth3) = 3
alg_order(alg::McAte3) = 3
alg_order(alg::McAte4) = 4
alg_order(alg::CandyRoz4) = 4
alg_order(alg::CalvoSanz4) = 4
alg_order(alg::McAte42) = 4
alg_order(alg::McAte5) = 5
alg_order(alg::Yoshida6) = 6
alg_order(alg::KahanLi6) = 6
alg_order(alg::McAte8) = 8
alg_order(alg::KahanLi8) = 8
alg_order(alg::SofSpa10) = 10

alg_order(alg::IRKN3) = 3
alg_order(alg::Nystrom4) = 4
alg_order(alg::Nystrom4VelocityIndependent) = 4
alg_order(alg::IRKN4) = 4
alg_order(alg::Nystrom5VelocityIndependent) = 5
alg_order(alg::DPRKN6) = 6
alg_order(alg::DPRKN8) = 8
alg_order(alg::DPRKN12) = 12
alg_order(alg::ERKN4) = 4
alg_order(alg::ERKN5) = 5

alg_order(alg::Midpoint) = 2
alg_order(alg::GenericIIF1) = 1
alg_order(alg::GenericIIF2) = 2
alg_order(alg::CarpenterKennedy2N54) = 4
alg_order(alg::SSPRK22) = 2
alg_order(alg::SSPRK33) = 3
alg_order(alg::SSPRK53) = 3
alg_order(alg::SSPRK63) = 3
alg_order(alg::SSPRK73) = 3
alg_order(alg::SSPRK83) = 3
alg_order(alg::SSPRK432) = 3
alg_order(alg::SSPRK932) = 3
alg_order(alg::SSPRK54) = 4
alg_order(alg::SSPRK104) = 4
alg_order(alg::RK4) = 4
alg_order(alg::ExplicitRK) = alg.tableau.order

alg_order(alg::BS3) = 3
alg_order(alg::BS5) = 5
alg_order(alg::OwrenZen3) = 3
alg_order(alg::OwrenZen4) = 4
alg_order(alg::OwrenZen5) = 5

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
alg_order(alg::GenericImplicitEuler) = 1
alg_order(alg::GenericTrapezoid) = 2
alg_order(alg::ImplicitEuler) = 1
alg_order(alg::LinearImplicitEuler) = 1
alg_order(alg::MidpointSplitting) = 2
alg_order(alg::Trapezoid) = 2
alg_order(alg::ImplicitMidpoint) = 2
alg_order(alg::TRBDF2) = 2
alg_order(alg::SSPSDIRK2) = 2
alg_order(alg::SDIRK2) = 2
alg_order(alg::Kvaerno3) = 3
alg_order(alg::Kvaerno4) = 4
alg_order(alg::Kvaerno5) = 5
alg_order(alg::KenCarp3) = 3
alg_order(alg::KenCarp4) = 4
alg_order(alg::KenCarp5) = 5
alg_order(alg::Cash4) = 4
alg_order(alg::Hairer4) = 4
alg_order(alg::Hairer42) = 4
alg_order(alg::Feagin10) = 10
alg_order(alg::Feagin12) = 12
alg_order(alg::Feagin14) = 14

alg_order(alg::Rosenbrock23) = 2
alg_order(alg::Rosenbrock32) = 3
alg_order(alg::ROS3P) = 3
alg_order(alg::Rodas3) = 3
alg_order(alg::RosShamp4) = 4
alg_order(alg::Veldd4) = 4
alg_order(alg::Velds4) = 4
alg_order(alg::GRK4T) = 4
alg_order(alg::GRK4A) = 4
alg_order(alg::Ros4LStab) = 4
alg_order(alg::Rodas4) = 4
alg_order(alg::Rodas42) = 4
alg_order(alg::Rodas4P) = 4
alg_order(alg::Rodas5) = 5

alg_order(alg::AB3) = 3
alg_order(alg::AB4) = 4
alg_order(alg::AB5) = 5
alg_order(alg::ABM32) = 3
alg_order(alg::ABM43) = 4
alg_order(alg::ABM54) = 5

alg_order(alg::VCAB3) = 3
alg_order(alg::VCAB4) = 4
alg_order(alg::VCAB5) = 5
alg_order(alg::VCABM3) = 3
alg_order(alg::VCABM4) = 4
alg_order(alg::VCABM5) = 5

alg_order(alg::VCABM) = 1  #dummy value

alg_order(alg::IMEXEuler) = 1
alg_order(alg::CNAB2) = 2
alg_order(alg::CNLF2) = 2

alg_order(alg::AN5) = 5
alg_order(alg::JVODE) = 1  #dummy value

alg_order(alg::ABDF2) = 2
alg_order(alg::QNDF1) = 1
alg_order(alg::QNDF2) = 2

alg_order(alg::QNDF) = 1 #dummy value

alg_order(alg::ROCK2) = 2

alg_maximum_order(alg) = alg_order(alg)
alg_maximum_order(alg::CompositeAlgorithm) = maximum(alg_order(x) for x in alg.algs)

alg_adaptive_order(alg::ExplicitRK) = alg.tableau.adaptiveorder
alg_adaptive_order(alg::OrdinaryDiffEqAlgorithm) = alg_order(alg)-1
alg_adaptive_order(alg::DP8) = 6
alg_adaptive_order(alg::Feagin10) = 8
alg_adaptive_order(alg::Feagin12) = 10
alg_adaptive_order(alg::Feagin14) = 12

alg_adaptive_order(alg::Rosenbrock23) = 3
alg_adaptive_order(alg::Rosenbrock32) = 2

alg_adaptive_order(alg::GenericImplicitEuler) = 0
alg_adaptive_order(alg::GenericTrapezoid) = 1
alg_adaptive_order(alg::ImplicitEuler) = 0
alg_adaptive_order(alg::Trapezoid) = 1
# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better
alg_adaptive_order(alg::ImplicitMidpoint) = 1
# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better

alg_adaptive_order(alg::AN5) = 5

beta2_default(alg::OrdinaryDiffEqAlgorithm) = 2//(5alg_order(alg))
beta2_default(alg::FunctionMap) = 0
beta2_default(alg::DP8) = 0//1
beta2_default(alg::DP5) = 4//100
beta2_default(alg::DP5Threaded) = 4//100

beta1_default(alg::OrdinaryDiffEqAlgorithm,beta2) = 7//(10alg_order(alg))
beta1_default(alg::FunctionMap,beta2) = 0
beta1_default(alg::DP8,beta2) = typeof(beta2)(1//alg_order(alg)) - beta2/5
beta1_default(alg::DP5,beta2) = typeof(beta2)(1//alg_order(alg)) - 3beta2/4
beta1_default(alg::DP5Threaded,beta2) = typeof(beta2)(1//alg_order(alg)) - 3beta2/4

gamma_default(alg::OrdinaryDiffEqAlgorithm) = 9//10

qsteady_min_default(alg::OrdinaryDiffEqAlgorithm) = 1
qsteady_max_default(alg::OrdinaryDiffEqAlgorithm) = 1
qsteady_max_default(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = 6//5
# But don't re-use Jacobian if not adaptive: too risky and cannot pull back
qsteady_max_default(alg::OrdinaryDiffEqImplicitAlgorithm) = 1//1
qsteady_max_default(alg::AN5) = 3//2
qsteady_max_default(alg::JVODE) = 3//2
qsteady_max_default(alg::QNDF1) = 2//1
qsteady_max_default(alg::QNDF2) = 2//1
qsteady_max_default(alg::QNDF) = 2//1

FunctionMap_scale_by_time(alg::FunctionMap{scale_by_time}) where {scale_by_time} = scale_by_time

# SSP coefficients
"""
    ssp_coefficient(alg)

Return the SSP coefficient of the ODE algorithm `alg`. If one time step of size
`dt` with `alg` can be written as a convex combination of explicit Euler steps
with step sizes `cᵢ * dt`, the SSP coefficient is the minimal value of `1/cᵢ`.

# Examples
```julia-repl
julia> ssp_coefficient(SSPRK104())
6
```
"""
ssp_coefficient(alg) = error("$alg is not a strong stability preserving method.")
ssp_coefficient(alg::Euler) = 1
ssp_coefficient(alg::SSPRK22) = 1
ssp_coefficient(alg::SSPRK33) = 1
ssp_coefficient(alg::SSPRK53) = 2.65
ssp_coefficient(alg::SSPRK63) = 3.518
ssp_coefficient(alg::SSPRK73) = 4.2879
ssp_coefficient(alg::SSPRK83) = 5.107
ssp_coefficient(alg::SSPRK432) = 2
ssp_coefficient(alg::SSPRK932) = 6
ssp_coefficient(alg::SSPRK54) = 1.508
ssp_coefficient(alg::SSPRK104) = 6

# We shouldn't do this probably.
#ssp_coefficient(alg::ImplicitEuler) = Inf
ssp_coefficient(alg::SSPSDIRK2) = 4

# stability regions
alg_stability_size(alg::DP5) = 3.3066
alg_stability_size(alg::Tsit5) = 3.5068
alg_stability_size(alg::Vern6) = 4.8553
alg_stability_size(alg::Vern7) = 4.6400
alg_stability_size(alg::Vern8) = 5.8641
alg_stability_size(alg::Vern9) = 4.4762

alg_can_repeat_jac(alg::OrdinaryDiffEqAlgorithm) = false
alg_can_repeat_jac(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = true

function unwrap_alg(integrator, is_stiff)
  alg = integrator.alg
  iscomp = typeof(alg) <: CompositeAlgorithm
  if !iscomp
    return alg
  elseif typeof(alg.choice_function) <: AutoSwitch
    num = is_stiff ? 2 : 1
    return alg.algs[num]
  else
    return alg.algs[integrator.cache.current]
  end
end
