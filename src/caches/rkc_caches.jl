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
