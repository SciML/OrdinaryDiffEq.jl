mutable struct ROCK2ConstantCache{T,T2,zType} <: OrdinaryDiffEqConstantCache
  ms::SVector{46, Int}
  fp1::SVector{46, T}
  fp2::SVector{46, T}
  recf::Vector{T2}
  zprev::zType
  mdegprev::Int
  mdeg::Int
  recind::Int
end
@cache struct ROCK2Cache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  gprev::uType
  gprev2::uType
  tmp::uType
  atmp::uNoUnitsType
  fsalfirst::rateType
  k::rateType
  k2::rateType
  constantcache::ROCK2ConstantCache
end

function alg_cache(alg::ROCK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  constantcache = ROCK2ConstantCache(uEltypeNoUnits, uEltypeNoUnits, u) # WIP: not sure about what type to use in here
  gprev = similar(u)
  gprev2 = similar(u)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k2 = zero(rate_prototype)
  ROCK2Cache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, k2, constantcache)
end

function alg_cache(alg::ROCK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ROCK2ConstantCache(uEltypeNoUnits, uEltypeNoUnits, u) # WIP: not sure about what type to use in here
end

mutable struct ROCK4ConstantCache{T,T2,T3,T4,zType} <: OrdinaryDiffEqConstantCache
  ms::SVector{50, Int}
  fpa::Vector{T}
  fpb::Vector{T2}
  fpbe::Vector{T3}
  recf::Vector{T4}
  zprev::zType
  mdegprev::Int
  mdeg::Int
  recind::Int
end
@cache struct ROCK4Cache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  gprev::uType
  gprev2::uType
  gprev3::uType
  gprev4::uType
  gprev5::uType
  tmp::uType
  atmp::uNoUnitsType
  fsalfirst::rateType
  k::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  constantcache::ROCK4ConstantCache
end

function alg_cache(alg::ROCK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  constantcache = ROCK4ConstantCache(real(uEltypeNoUnits), real(uEltypeNoUnits), u)
  gprev = similar(u)
  gprev2 = similar(u)
  gprev3 = similar(u)
  gprev4 = similar(u)
  gprev5 = similar(u)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  ROCK4Cache(u, uprev, gprev, gprev2, gprev3, gprev4, gprev5, tmp, atmp, fsalfirst, k, k2, k3, k4, k5, constantcache)
end

function alg_cache(alg::ROCK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ROCK4ConstantCache(real(uEltypeNoUnits), real(uEltypeNoUnits),u)
end

mutable struct RKCConstantCache{zType} <: OrdinaryDiffEqConstantCache
  #to match the types to call maxeig!
  zprev::zType
end
@cache struct RKCCache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  gprev::uType
  gprev2::uType
  tmp::uType
  atmp::uNoUnitsType
  fsalfirst::rateType
  k::rateType
  constantcache::RKCConstantCache
end

function alg_cache(alg::RKC,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  constantcache = RKCConstantCache(u)
  gprev = similar(u)
  gprev2 = similar(u)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  RKCCache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::RKC,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  RKCConstantCache(u)
end

@cache mutable struct IRKCConstantCache{uType,rateType,F,N} <: OrdinaryDiffEqConstantCache
  minm::Int64
  zprev::uType
  k2::rateType
  uf::F
  nlsolve::N
  du₁::rateType
  du₂::rateType
end

@cache mutable struct IRKCCache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  gprev::uType
  gprev2::uType
  fsalfirst::rateType
  k::rateType
  du1::rateType
  f1ⱼ₋₁::rateType
  f1ⱼ₋₂::rateType
  f2ⱼ₋₁::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
  utilde::uType
  atmp::uNoUnitsType
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  du₁::rateType
  du₂::rateType
  constantcache::IRKCConstantCache
end

function alg_cache(alg::IRKC,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//1, 1//1
  @oopnlsolve
  zprev = u
  k2  = rate_prototype
  du₁ = rate_prototype; du₂ = rate_prototype
  IRKCConstantCache(50,zprev,k2,uf,nlsolve,du₁,du₂)
end

function alg_cache(alg::IRKC,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//1, 1//1
  @iipnlsolve

  gprev = similar(u)
  gprev2 = similar(u)
  tmp = similar(u)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  fsalfirst = zero(rate_prototype)
  k  = zero(rate_prototype)
  k2 = zero(rate_prototype)
  zprev = similar(u)
  f1ⱼ₋₁ = zero(rate_prototype)
  f1ⱼ₋₂ = zero(rate_prototype)
  f2ⱼ₋₁ = zero(rate_prototype)
  du₁ = zero(rate_prototype)
  du₂ = zero(rate_prototype)
  constantcache = IRKCConstantCache(50,zprev,k2,uf,nlsolve,du₁,du₂)
  IRKCCache(u,uprev,gprev,gprev2,fsalfirst,k,du1,f1ⱼ₋₁,f1ⱼ₋₂,f2ⱼ₋₁,z,dz,b,tmp,utilde,atmp,J,W,uf,jac_config,linsolve,nlsolve,du₁,du₂,constantcache)
end
