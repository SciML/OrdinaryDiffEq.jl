# IMEXEuler

mutable struct IMEXEulerConstantCache{F,uToltype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

mutable struct IMEXEulerCache{uType,rateType,uNoUnitsType,J,W,UF,JC,uToltype,F} <: OrdinaryDiffEqMutableCache
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
  if DiffEqBase.has_jac(f)
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  z = similar(u,axes(u))
  dz = similar(u,axes(u)); tmp = similar(u,axes(u)); b = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  du₁ = zero(rate_prototype)
  du1 = zero(rate_prototype)
  atmp = similar(u,uEltypeNoUnits,axes(u))

  uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  uToltype = real(uBottomEltypeNoUnits)
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
  IMEXEulerCache(u,uprev,uprev2,du1,fsalfirst,k,du₁,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end
