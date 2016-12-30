abstract OrdinaryDiffEqCache <: DECache
immutable ODEEmptyCache <: OrdinaryDiffEqCache end

alg_cache{F}(alg::OrdinaryDiffEqAlgorithm,prob,callback::F) = ODEEmptyCache()

immutable EulerCache{uType,rateType} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  k::rateType
end

function alg_cache{uType<:AbstractArray}(alg::Euler,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  EulerCache(u,uprev,similar(rate_prototype))
end

alg_cache{uType}(alg::Euler,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable MidpointCache{uType,rateType} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  k::rateType
  kprev::rateType
  du::rateType
  utilde::uType
end

function alg_cache{uType<:AbstractArray}(alg::Midpoint,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  utilde = similar(u)
  k = similar(rate_prototype)
  du = similar(rate_prototype)
  MidpointCache(u,uprev,kprev,k,du,utilde)
end

alg_cache{uType}(alg::Midpoint,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable RK4Cache{uType,rateType} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::rateType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  tmp::uType
end

function alg_cache{uType<:AbstractArray}(alg::RK4,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k₁ = similar(rate_prototype)
  k₂ = similar(rate_prototype)
  k₃ = similar(rate_prototype)
  k₄ = similar(rate_prototype)
  tmp = similar(u)
  RK4Cache(u,uprev,kprev,k₁,k₂,k₃,k₄,tmp)
end

alg_cache{uType}(alg::RK4,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable BS3Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::rateType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

function alg_cache{uType<:AbstractArray}(alg::BS3,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS3Cache(u,uprev,kprev,k1,k2,k3,k4,utilde,tmp,atmp)
end

alg_cache{uType}(alg::BS3,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable BS5Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::Vector{rateType}
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  utilde::uType
  uhat::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmptilde::uEltypeNoUnits
end

function alg_cache{uType<:AbstractArray}(alg::BS5,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  k8 = similar(rate_prototype)
  utilde = similar(u); uhat = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits)
  BS5Cache(u,uprev,kprev,k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde)
end

alg_cache{uType}(alg::BS5,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable Tsit5Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::Vector{rateType}
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

function alg_cache{uType<:AbstractArray}(alg::Tsit5,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits)
  Tsit5Cache(u,uprev,kprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp)
end

alg_cache{uType}(alg::Tsit5,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable ExplicitRKCache{uType,rateType,uEltypeNoUnits,ksEltype} <: OrdinaryDiffEqCache
  u::uType
  tmp::uType
  utilde::rateType
  uEEst::rateType
  atmp::uEltypeNoUnits
  uprev::uType
  kprev::ksEltype
  utmp::uType
  kk::Vector{ksEltype}
end

function alg_cache{uType<:AbstractArray}(alg::ExplicitRK,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  kk = Vector{typeof(rate_prototype)}(0)
  for i = 1:tableau.stages
    push!(kk,similar(rate_prototype))
  end
  utilde = similar(rate_prototype)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  utmp = zeros(u)
  uEEst = similar(rate_prototype)
  ExplicitRKCache(u,tmp,utilde,uEEst,atmp,uprev,kprev,utmp,kk)
end

alg_cache{uType}(alg::ExplicitRK,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable Feagin10Cache{uType,uEltypeNoUnits,rateType} <: OrdinaryDiffEqCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  utmp::uType
  uprev::uType
  kprev::rateType
end

function alg_cache{uType<:AbstractArray}(alg::Feagin10,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  utmp = similar(u);

  Feagin10Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,utmp,uprev,kprev)
end

alg_cache{uType}(alg::Feagin10,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()

immutable Feagin12Cache{uType,uEltypeNoUnits,rateType} <: OrdinaryDiffEqCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  k18::rateType
  k19::rateType
  k20::rateType
  k21::rateType
  k22::rateType
  k23::rateType
  k24::rateType
  k25::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  utmp::uType
  uprev::uType
  kprev::rateType
end

function alg_cache{uType<:AbstractArray}(alg::Feagin12,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype); k18 = similar(rate_prototype)
  k19 = similar(rate_prototype); k20 = similar(rate_prototype); k21 = similar(rate_prototype); k22 = similar(rate_prototype)
  k23 = similar(rate_prototype); k24 = similar(rate_prototype); k25 = similar(rate_prototype)
  utmp = similar(u);
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)

  Feagin12Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,utmp,uprev,kprev)
end

alg_cache{uType}(alg::Feagin12,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()


immutable Feagin14Cache{uType,uEltypeNoUnits,rateType} <: OrdinaryDiffEqCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  k18::rateType
  k19::rateType
  k20::rateType
  k21::rateType
  k22::rateType
  k23::rateType
  k24::rateType
  k25::rateType
  k26::rateType
  k27::rateType
  k28::rateType
  k29::rateType
  k30::rateType
  k31::rateType
  k32::rateType
  k33::rateType
  k34::rateType
  k35::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  utmp::uType
  uprev::uType
  kprev::rateType
end

function alg_cache{uType<:AbstractArray}(alg::Feagin14,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype); k18 = similar(rate_prototype)
  k19 = similar(rate_prototype); k20 = similar(rate_prototype); k21 = similar(rate_prototype); k22 = similar(rate_prototype)
  k23 = similar(rate_prototype); k24 = similar(rate_prototype); k25 = similar(rate_prototype)
  k26 = similar(rate_prototype); k27 = similar(rate_prototype); k28 = similar(rate_prototype)
  k29 = similar(rate_prototype); k30 = similar(rate_prototype); k31 = similar(rate_prototype); k32 = similar(rate_prototype)
  k33 = similar(rate_prototype); k34 = similar(rate_prototype); k35 = similar(rate_prototype)
  utmp = similar(u);
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)

  Feagin14Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
                k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,
                k31,k32,k33,k34,k35,tmp,atmp,utmp,uprev,kprev)
end

alg_cache{uType}(alg::Feagin14,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()
