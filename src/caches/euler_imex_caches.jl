# IMEXEuler

mutable struct IMEXEulerConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

mutable struct IMEXEulerCache{uType,rateType,uNoUnitsType,J,W,UF,JC,N,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  du₁::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::J
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

u_cache(c::IMEXEulerCache)    = ()
du_cache(c::IMEXEulerCache)   = ()

function alg_cache(alg::IMEXEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  uf != nothing && ( uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p) )
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))
  IMEXEulerConstantCache(uf,nlsolve)
end

function alg_cache(alg::IMEXEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  du₁ = zero(rate_prototype)
  atmp = similar(u,uEltypeNoUnits,axes(u))
  if uf != nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
    linsolve = alg.linsolve(Val{:init},uf,u)
    jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  end
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))
  IMEXEulerCache(u,uprev,uprev2,du1,fsalfirst,k,du₁,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve)
end
