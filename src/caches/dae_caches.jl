@cache mutable struct DImplicitEulerCache{uType,uNoUnitsType,N} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  tmp::uType
  atmp::uNoUnitsType
  nlsolver::N
end

mutable struct DImplicitEulerConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::DImplicitEuler,du,u,res_prototype,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  DImplicitEulerConstantCache()
end

function alg_cache(alg::DImplicitEuler,du,u,res_prototype,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  γ, c = 1, 1
  α = 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,res_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,α,Val(true))

  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)

  DImplicitEulerCache(u,uprev,uprev2,tmp,atmp,nlsolver)
end
