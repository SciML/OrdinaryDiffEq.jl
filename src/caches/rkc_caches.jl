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
struct ROCK2Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  gprev::uType
  gprev2::uType
  tmp::uType
  atmp::uEltypeNoUnits
  fsalfirst::rateType
  k::rateType
  k2::rateType
  constantcache::ROCK2ConstantCache
end
u_cache(c::ROCK2Cache)    = (c.atmp,)
du_cache(c::ROCK2Cache)   = (c.fsalfirst,c.k,c.k2)

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
