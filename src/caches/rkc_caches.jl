struct ROCK2Cache{} <: OrdinaryDiffEqMutableCache # WIP
  constantcache::ROCK2ConstantCache
end

function alg_cache(alg::ROCK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  ROCK2Cache()
end

function alg_cache(alg::ROCK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  ROCK2ConstantCache(uEltypeNoUnits, uEltypeNoUnits, u) # WIP: not sure about what type to use in here
end
