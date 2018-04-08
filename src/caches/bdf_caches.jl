mutable struct ABDF2ConstantCache{F,uToltype,dtType,rate_prototype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  eulercache::ImplicitEulerConstantCache
  dtₙ₊₁::dtType
  dtₙ::dtType
  uₙ₋₁::rate_prototype
end

function alg_cache(alg::ABDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
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

  eulercache = ImplicitEulerConstantCache(uf,ηold,κ,tol,100000)

  dtₙ₊₁ = one(dt)
  dtₙ = one(dt)
  uₙ₋₁ = zero(u)

  ABDF2ConstantCache(uf, ηold, κ, tol, 10000, eulercache,
                      dtₙ₊₁, dtₙ, uₙ₋₁)
end

#=
mutable struct ABDF2Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,F,dtType} <: OrdinaryDiffEqMutableCache
  du1::rateType
  fsalfirst::rateType
  f2::rateType
  zprev::uType
  zᵧ::uType
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
  dtprev::dtType
end

u_cache(c::ABDF2Cache)    = (c.zprev,c.zᵧ,c.z,c.dz)
du_cache(c::ABDF2Cache)   = (c.f2,c.fsalfirst)

function alg_cache(alg::ABDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  zprev = similar(u,indices(u))
  zᵧ = similar(u,indices(u)); z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zeros(rate_prototype)
  f2 = zeros(rate_prototype)
  tmp = similar(u); b = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))

  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
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

  dtprev = zero(dt)

  ABDF2Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              du1,fsalfirst,f2,zprev,zᵧ,z,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,dtprev)
end
=#
