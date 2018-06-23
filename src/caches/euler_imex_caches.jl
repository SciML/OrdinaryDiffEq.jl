# IMEXEuler

mutable struct IMEXEulerConstantCache{F,uToltype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

mutable struct IMEXEulerCache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,tType,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::J
  W::J
  uf::UF
  jac_config::JC
  linsolve::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

u_cache(c::IMEXEulerCache)    = ()
du_cache(c::IMEXEulerCache)   = ()

function alg_cache(alg::IMEXEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  ηold = one(uToltype)

  if alg.κ != nothing
    κ = uToltype(alg.κ)
  else
    κ = uToltype(1//100)
  end
  if alg.tol != nothing
    tol = uToltype(alg.tol)
  else
    tol = uToltype(min(0.03,first(reltol)^(0.5)))
  end

  IMEXEulerConstantCache(uf,ηold,κ,tol,10000)
end

function alg_cache(alg::IMEXEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  J = zeros(uEltypeNoUnits,length(u),length(u))
  W = similar(J)
  z = similar(u,indices(u))
  dz = similar(u,indices(u))
  tmp = similar(u); b = similar(u,indices(u));
  atmp = similar(u,uEltypeNoUnits,indices(u))
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du₁ = zeros(rate_prototype)

  if typeof(f) <: SplitFunction
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end

  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)

  uToltype = real(uBottomEltypeNoUnits)
  if alg.κ != nothing
    κ = uToltype(alg.κ)
  else
    κ = uToltype(1//100)
  end
  if alg.tol != nothing
    tol = uToltype(alg.tol)
  else
    tol = uToltype(min(0.03,first(reltol)^(0.5)))
  end

  uprev3 = similar(u)
  tprev2 = t

  ηold = one(uToltype)

  IMEXEulerCache(u,uprev,uprev2,fsalfirst,k,k1,k2,du₁,du1,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000,uprev3,tprev2)
end
