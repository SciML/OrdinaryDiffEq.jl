# SBDF2

mutable struct SBDF2ConstantCache{rateType,F,uToltype,uType} <: OrdinaryDiffEqConstantCache
  k2::rateType
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  uprev2::uType
  k₁::rateType
  k₂::rateType
end

mutable struct SBDF2Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,uToltype,tType,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k::rateType
  du1::rateType
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
  uprev2::uType
  k₁::rateType
  k₂::rateType
  du₁::rateType
end

u_cache(c::SBDF2Cache)    = ()
du_cache(c::SBDF2Cache)   = ()

function alg_cache(alg::SBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  k2 = rate_prototype
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  ηold = one(uToltype)
  uprev2 = u
  k₁ = rate_prototype
  k₂ = rate_prototype

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

  SBDF2ConstantCache(k2,uf,ηold,κ,tol,10000,uprev2,k₁,k₂)
end

function alg_cache(alg::SBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  z = similar(u,axes(u))
  dz = similar(u,axes(u))
  tmp = similar(u); b = similar(u,axes(u));
  atmp = similar(u,uEltypeNoUnits,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  du1 = zero(rate_prototype)

  uprev2 = similar(u)
  k₁ = zero(rate_prototype)
  k₂ = zero(rate_prototype)
  du₁ = zero(rate_prototype)

  uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)

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

  ηold = one(uToltype)

  SBDF2Cache(u,uprev,fsalfirst,k,du1,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000,uprev2,k₁,k₂,du₁)
end