@cache struct Rice3ConstantCache{TabType} <: OrdinaryDiffEqConstantCache
  tab::TabType
end

@cache struct Rice3Cache{uType,fsalType,rateType,rateTypeFast,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::fsalType
  k2::rateType
  k3::rateType
  k4::rateType
  h2::rateTypeFast
  h3::rateTypeFast
  h4::rateTypeFast
  d1::rateTypeFast
  d2::rateTypeFast
  d3::rateTypeFast
  utilde::uType
  tmp::uType
  tmp2::uType
  atmp::uNoUnitsType
  tab::TabType
end

function alg_cache(alg::Rice3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = alg_cache(alg.base_method,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val(false))
  k1 = zero(rate_prototype.x[1])
  k2 = zero(rate_prototype.x[1])
  k3 = zero(rate_prototype.x[1])
  k4 = zero(rate_prototype.x[1])
  h1 = zero(rate_prototype.x[2])
  h2 = zero(rate_prototype.x[2])
  h3 = zero(rate_prototype.x[2])
  h4 = zero(rate_prototype.x[2])
  d1 = zero(rate_prototype.x[2])
  d2 = zero(rate_prototype.x[2])
  d3 = zero(rate_prototype.x[2])
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  tmp2 = similar(u)
  Rice3Cache(u,uprev,ArrayPartition(k1,h1),k2,k3,k4,h2,h3,h4,d1,d2,d3,utilde,tmp,tmp2,atmp,tab)
end

function alg_cache(alg::Rice3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = alg_cache(alg.base_method,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,Val(false))
  return Rice3ConstantCache(tab)
end
