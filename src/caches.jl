abstract OrdinaryDiffEqCache <: DECache
immutable ODEEmptyCache <: OrdinaryDiffEqCache end
immutable ODEChunkCache{CS} <: OrdinaryDiffEqCache end

alg_cache{F}(alg::OrdinaryDiffEqAlgorithm,prob,callback::F) = ODEEmptyCache()

immutable EulerCache{uType,rateType} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  k::rateType
  fsalfirst::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Euler,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  EulerCache(u,uprev,similar(rate_prototype),similar(rate_prototype))
end

alg_cache{uType}(alg::Euler,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable MidpointCache{uType,rateType} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  k::rateType
  kprev::rateType
  du::rateType
  utilde::uType
  fsalfirst::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Midpoint,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  utilde = similar(u)
  k = similar(rate_prototype)
  du = similar(rate_prototype)
  fsalfirst = similar(rate_prototype)
  MidpointCache(u,uprev,kprev,k,du,utilde,fsalfirst)
end

alg_cache{uType}(alg::Midpoint,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable RK4Cache{uType,rateType} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::rateType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::RK4,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k₁ = similar(rate_prototype)
  k₂ = similar(rate_prototype)
  k₃ = similar(rate_prototype)
  k₄ = similar(rate_prototype)
  k = similar(rate_prototype)
  tmp = similar(u)
  RK4Cache(u,uprev,kprev,k₁,k₂,k₃,k₄,k,tmp)
end

alg_cache{uType}(alg::RK4,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

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

Base.@pure function alg_cache{uType<:AbstractArray}(alg::BS3,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS3Cache(u,uprev,kprev,k1,k2,k3,k4,utilde,tmp,atmp)
end

alg_cache{uType}(alg::BS3,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable BS5Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
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

Base.@pure function alg_cache{uType<:AbstractArray}(alg::BS5,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
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
  BS5Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde)
end

alg_cache{uType}(alg::BS5,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable Tsit5Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
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

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Tsit5,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
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
  Tsit5Cache(u,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp)
end

alg_cache{uType}(alg::Tsit5,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable DP5Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  update::rateType
  bspl::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::DP5,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  dense_tmp3 = similar(rate_prototype)
  dense_tmp4 = similar(rate_prototype)
  update = similar(rate_prototype)
  bspl = similar(rate_prototype)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  DP5Cache(u,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp)
end

alg_cache{uType}(alg::DP5,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

alg_cache(alg::DP5Threaded,u,rate_prototype,uEltypeNoUnits,uprev,kprev) = alg_cache(DP5(),u,rate_prototype,uEltypeNoUnits,uprev,kprev)

immutable Vern6Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
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
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Vern6,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype);
  k8 = similar(rate_prototype); k9 = similar(rate_prototype);
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  Vern6Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,atmp)
end

alg_cache{uType}(alg::Vern6,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable Vern7Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
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
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Vern7,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype);
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); utilde = similar(u); update = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  Vern7Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,update,tmp,atmp)
end

alg_cache{uType}(alg::Vern7,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()


immutable Vern8Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
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
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Vern8,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype);
  k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype);
  k8 = similar(rate_prototype); tmp = similar(u)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype);
  k12 = similar(rate_prototype); k13 = similar(rate_prototype)
  utilde = similar(u); update = similar(u);
  atmp = similar(u,uEltypeNoUnits)
  Vern8Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp)
end

alg_cache{uType}(alg::Vern8,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable Vern9Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
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
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Vern9,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype);k3 = similar(rate_prototype);
  k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype);k7 = similar(rate_prototype);
  k8 = similar(rate_prototype);
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype);
  k12 = similar(rate_prototype); update = similar(u)
  k13 = similar(rate_prototype); k14 = similar(rate_prototype); k15 = similar(rate_prototype);
  k16 =similar(rate_prototype);
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  Vern9Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,update,tmp,atmp)
end

alg_cache{uType}(alg::Vern9,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable TanYam7Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::rateType
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
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::TanYam7,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype) ; k3 = similar(rate_prototype); k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6 = similar(rate_prototype) ; k7 = similar(rate_prototype); k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10= similar(rate_prototype) ;
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)
  TanYam7Cache(u,uprev,kprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k)
end

alg_cache{uType}(alg::TanYam7,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()


immutable DP8Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
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
  kupdate::rateType
  update::uType
  udiff::rateType
  bspl::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  dense_tmp5::rateType
  dense_tmp6::rateType
  dense_tmp7::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmp2::uEltypeNoUnits
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::DP8,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2  = similar(rate_prototype); k3  = similar(rate_prototype);  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6  = similar(rate_prototype); k7  = similar(rate_prototype);  k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype); k12 = similar(rate_prototype)
  kupdate = similar(rate_prototype); utilde = similar(u);
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); atmp2 = similar(u,uEltypeNoUnits); update = similar(u)
  k13 = similar(rate_prototype)
  k14 = similar(rate_prototype)
  k15 = similar(rate_prototype)
  k16 = similar(rate_prototype)
  udiff = similar(rate_prototype)
  bspl = similar(rate_prototype)
  # dense_tmp1 = udiff
  # dense_tmp2 = bspl
  dense_tmp3 = similar(rate_prototype)
  dense_tmp4 = similar(rate_prototype)
  dense_tmp5 = similar(rate_prototype)
  dense_tmp6 = similar(rate_prototype)
  dense_tmp7 = similar(rate_prototype)
  DP8Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,kupdate,
           update,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,
           utilde,tmp,atmp,atmp2)
end

alg_cache{uType}(alg::DP8,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable TsitPap8Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCache
  u::uType
  uprev::uType
  kprev::rateType
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
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::TsitPap8,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype); k12 = similar(rate_prototype)
  k13 = similar(rate_prototype); update = similar(u); utilde = similar(u); k = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  TsitPap8Cache(u,uprev,kprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,k)
end

alg_cache{uType}(alg::TsitPap8,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable ExplicitRKCache{uType,rateType,uEltypeNoUnits,ksEltype} <: OrdinaryDiffEqCache
  u::uType
  tmp::uType
  utilde::rateType
  uEEst::rateType
  atmp::uEltypeNoUnits
  uprev::uType
  kprev::ksEltype
  kk::Vector{ksEltype}
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::ExplicitRK,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  kk = Vector{typeof(rate_prototype)}(0)
  for i = 1:alg.tableau.stages
    push!(kk,similar(rate_prototype))
  end
  utilde = similar(rate_prototype)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  uEEst = similar(rate_prototype)
  ExplicitRKCache(u,tmp,utilde,uEEst,atmp,uprev,kprev,kk)
end

alg_cache{uType}(alg::ExplicitRK,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

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
  uprev::uType
  kprev::rateType
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Feagin10,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)

  Feagin10Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,uprev,kprev,k)
end

alg_cache{uType}(alg::Feagin10,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

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
  uprev::uType
  kprev::rateType
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Feagin12,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype); k18 = similar(rate_prototype)
  k19 = similar(rate_prototype); k20 = similar(rate_prototype); k21 = similar(rate_prototype); k22 = similar(rate_prototype)
  k23 = similar(rate_prototype); k24 = similar(rate_prototype); k25 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)

  Feagin12Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,uprev,kprev,k)
end

alg_cache{uType}(alg::Feagin12,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()


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
  uprev::uType
  kprev::rateType
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Feagin14,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype); k18 = similar(rate_prototype)
  k19 = similar(rate_prototype); k20 = similar(rate_prototype); k21 = similar(rate_prototype); k22 = similar(rate_prototype)
  k23 = similar(rate_prototype); k24 = similar(rate_prototype); k25 = similar(rate_prototype)
  k26 = similar(rate_prototype); k27 = similar(rate_prototype); k28 = similar(rate_prototype)
  k29 = similar(rate_prototype); k30 = similar(rate_prototype); k31 = similar(rate_prototype); k32 = similar(rate_prototype)
  k33 = similar(rate_prototype); k34 = similar(rate_prototype); k35 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)

  Feagin14Cache(u,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
                k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,
                k31,k32,k33,k34,k35,tmp,atmp,uprev,kprev,k)
end

alg_cache{uType}(alg::Feagin14,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()


immutable LowOrderRosenbrockCache{uType,rateType,vecuType,JType} <: OrdinaryDiffEqCache
  u::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uType
  J::JType
  W::JType
  tmp2::uType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Rosenbrock23,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  k₁ = similar(rate_prototype)
  k₂ = similar(rate_prototype)
  k₃ = similar(rate_prototype)
  du1 = similar(rate_prototype)
  du2 = similar(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = similar(rate_prototype)
  vectmp = similar(vec(u))
  vectmp2 = similar(vec(u))
  vectmp3 = similar(vec(u))
  fsalfirst = similar(rate_prototype)
  fsallast = similar(rate_prototype)
  dT = similar(u)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp2 = similar(u)

  LowOrderRosenbrockCache(u,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp2)
end

alg_cache{uType<:AbstractArray}(alg::Rosenbrock32,u::uType,
                                rate_prototype,uEltypeNoUnits,
                                uprev,kprev) =
                                alg_cache(Rosenbrock23(),u,rate_prototype,uEltypeNoUnits,uprev,kprev)

alg_cache{uType}(alg::Rosenbrock23,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()
alg_cache{uType}(alg::Rosenbrock32,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev) = ODEEmptyCache()

immutable ImplicitEulerCache{uType,DiffCacheType,rateType,CS} <: OrdinaryDiffEqCache
  u::uType
  dual_cache::DiffCacheType
  u_old::uType
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::ImplicitEuler,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  u_old = similar(u); k = similar(rate_prototype)
  if get_chunksize(alg) != 0
    dual_cache = DiffCache(u,alg)
    ImplicitEulerCache{typeof(u),typeof(dual_cache),typeof(k),get_chunksize(alg)}(u,dual_cache,u_old,k)
  else
    dual_cache = DiffCache(u)
    ImplicitEulerCache{typeof(u),typeof(dual_cache),typeof(k),ForwardDiff.pickchunksize(length(u))}(u,dual_cache,u_old,k)
  end
end

Base.@pure function alg_cache{uType}(alg::ImplicitEuler,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  if get_chunksize(alg) != 0
    return ODEChunkCache{get_chunksize(alg)}()
  else
    return ODEChunkCache{ForwardDiff.pickchunksize(length(u))}()
  end
end

immutable TrapezoidCache{uType,DiffCacheType,rateType,CS} <: OrdinaryDiffEqCache
  u::uType
  dual_cache::DiffCacheType
  dual_cache2::DiffCacheType
  u_old::uType
  k::rateType
end

Base.@pure function alg_cache{uType<:AbstractArray}(alg::Trapezoid,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  u_old = similar(u);k = similar(rate_prototype)
  if get_chunksize(alg) != 0
    dual_cache = DiffCache(u,alg)
    dual_cache2 = DiffCache(u,alg)
    TrapezoidCache{typeof(u),typeof(dual_cache),typeof(k),get_chunksize(alg)}(u,dual_cache,dual_cache2,u_old,k)
  else
    dual_cache = DiffCache(u)
    dual_cache2 = DiffCache(u)
    TrapezoidCache{typeof(u),typeof(dual_cache),typeof(k),ForwardDiff.pickchunksize(length(u))}(u,dual_cache,dual_cache2,u_old,k)
  end
end

Base.@pure function alg_cache{uType}(alg::Trapezoid,u::uType,rate_prototype,uEltypeNoUnits,uprev,kprev)
  if get_chunksize(alg) != 0
    return ODEChunkCache{get_chunksize(alg)}()
  else
    return ODEChunkCache{ForwardDiff.pickchunksize(length(u))}()
  end
end

get_chunksize(cache::DECache) = error("This cache does not have a chunksize.")
Base.@pure get_chunksize{uType,DiffCacheType,rateType,CS}(cache::ImplicitEulerCache{uType,DiffCacheType,rateType,CS}) = CS
Base.@pure get_chunksize{uType,DiffCacheType,rateType,CS}(cache::TrapezoidCache{uType,DiffCacheType,rateType,CS}) = CS
Base.@pure get_chunksize{CS}(cache::ODEChunkCache{CS}) = CS
