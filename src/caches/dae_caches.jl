@cache mutable struct DImplicitEulerCache{uType,rateType,uNoUnitsType,N} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  atmp::uNoUnitsType
  nlsolver::N
end

function alg_cache(alg::DImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  atmp = similar(u,uEltypeNoUnits)

  DImplicitEulerCache(u,uprev,uprev2,fsalfirst,atmp,nlsolver)
end

mutable struct DImplicitEulerConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::DImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1, 1
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  DImplicitEulerConstantCache(nlsolver)
end
