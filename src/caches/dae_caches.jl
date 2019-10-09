@cache mutable struct DImplicitEulerCache{uType,uNoUnitsType,N} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  tmp::uType
  atmp::uNoUnitsType
  nlsolver::N
end

mutable struct DImplicitEulerConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::DImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  DImplicitEulerConstantCache()
end

function alg_cache(alg::DImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  γ, c = 1.0, 1.0
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))

  DImplicitEulerCache(u,uprev,uprev2,tmp,atmp,nlsolver)
end
