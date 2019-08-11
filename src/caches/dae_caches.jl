@cache mutable struct DImplicitEulerCache{uType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  tmp::uType
  atmp::uNoUnitsType
end

mutable struct DImplicitEulerConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::DImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  DImplicitEulerConstantCache()
end

function alg_cache(alg::DImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)

  DImplicitEulerCache(u,uprev,uprev2,tmp,atmp)
end
