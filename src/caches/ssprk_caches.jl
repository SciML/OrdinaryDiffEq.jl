immutable SSPRK22Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::SSPRK22Cache) = ()
du_cache(c::SSPRK22Cache) = (c.k,c.fsalfirst)

immutable SSPRK22ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK22Cache(u,uprev,k,tmp,fsalfirst)
end

alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK22ConstantCache()


immutable SSPRK33Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::SSPRK33Cache) = ()
du_cache(c::SSPRK33Cache) = (c.k,c.fsalfirst)

immutable SSPRK33ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK33Cache(u,uprev,k,tmp,fsalfirst)
end

alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK33ConstantCache()


immutable SSPRK432Cache{uType,rateType,uArrayType,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  utilde::uArrayType
  atmp::uEltypeNoUnits
end

u_cache(c::SSPRK432Cache) = (c.utilde,c.atmp)
du_cache(c::SSPRK432Cache) = (c.k,c.fsalfirst)

immutable SSPRK432ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  SSPRK432Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp)
end

alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK432ConstantCache()


immutable SSPRK104Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₄::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::SSPRK104Cache) = ()
du_cache(c::SSPRK104Cache) = (c.k,c.fsalfirst,c.k₄)

immutable SSPRK104ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK104Cache(u,uprev,k,k₄,tmp,fsalfirst)
end

alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK104ConstantCache()
