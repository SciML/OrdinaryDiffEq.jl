struct AB3Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  ralk2::rateType
  u1::uType
  u2::uType
  k::rateType
  tmp::uType
end

u_cache(c::AB3Cache) = ()
du_cache(c::AB3Cache) = (c.fsalfirst,c.k2,c.k3,c.ralk2,c.k)

struct AB3ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::AB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  ralk2 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  AB3Cache(u,uprev,k1,k2,k3,ralk2,u1,u2,k,tmp)
end

alg_cache(alg::AB3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = AB3ConstantCache()

struct ABM32Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  ralk2::rateType
  u2::uType
  k::rateType
  tmp::uType
end

u_cache(c::ABM32Cache) = ()
du_cache(c::ABM32Cache) = (c.fsalfirst,c.k2,c.k3,c.ralk2,c.k)

struct ABM32ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::ABM32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  ralk2 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  ABM32Cache(u,uprev,k1,k2,k3,ralk2,u2,k,tmp)
end

alg_cache(alg::ABM32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = ABM32ConstantCache()