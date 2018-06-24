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
struct ROCK2Cache{} <: OrdinaryDiffEqMutableCache # WIP
  constantcache::ROCK2ConstantCache
end

function alg_cache(alg::ROCK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  ROCK2Cache()
end

function alg_cache(alg::ROCK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ROCK2ConstantCache(uEltypeNoUnits, uEltypeNoUnits, u) # WIP: not sure about what type to use in here
end
