mutable struct ABDF2ConstantCache{F,uToltype,dtType,rate_prototype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  eulercache::ImplicitEulerConstantCache
  dtₙ₋₁::dtType
  fsalfirstprev::rate_prototype
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

  dtₙ₋₁ = one(dt)
  fsalfirstprev = rate_prototype

  ABDF2ConstantCache(uf, ηold, κ, tol, 10000, eulercache, dtₙ₋₁, fsalfirstprev)
end

mutable struct ABDF2Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,F,dtType} <: OrdinaryDiffEqMutableCache
  uₙ::uType
  uₙ₋₁::uType
  uₙ₋₂::uType
  du1::rateType
  fsalfirst::rateType
  fsalfirstprev::rateType
  k::rateType
  z::uType
  zₙ₋₁::uType
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
  eulercache::ImplicitEulerCache
  dtₙ₋₁::dtType
end

u_cache(c::ABDF2Cache)    = (c.z,c.dz)
du_cache(c::ABDF2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ABDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  zprev = similar(u,axes(u))
  zₙ₋₁ = similar(u,axes(u)); z = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  fsalfirstprev = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u))
  atmp = similar(u,uEltypeNoUnits,axes(u))

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

  eulercache = ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000)

  dtₙ₋₁ = one(dt)
  ABDF2Cache(u,uprev,uprev2,du1,fsalfirst,fsalfirstprev,k,z,zₙ₋₁,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,eulercache,dtₙ₋₁)
end

# QNDF1

mutable struct QNDF1ConstantCache{F,uToltype,coefType,coefType1,dtType,uType,tType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  eulercache::ImplicitEulerConstantCache
  D::coefType
  temp_D::coefType
  D2::coefType1
  R::coefType
  U::coefType
  uprev2::uType
  uprev3::uType
  tprev2::tType
  dtₙ₋₁::dtType
end

mutable struct QNDF1Cache{uType,rateType,coefType,uNoUnitsType,J,UF,JC,uToltype,F,dtType} <: OrdinaryDiffEqMutableCache
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  fsalfirstprev::rateType
  k::rateType
  z::uType
  zₙ₋₁::uType
  dz::uType
  b::uType
  D::coefType
  temp_D::coefType
  temp_u::uType
  D2::coefType
  R::coefType
  U::coefType
  tmp::uType
  atmp::uNoUnitsType
  utilde::uType
  J::J
  W::J
  uf::UF
  jac_config::JC
  linsolve::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  eulercache::ImplicitEulerCache
  dtₙ₋₁::dtType
end

u_cache(c::QNDF1Cache)    = ()
du_cache(c::QNDF1Cache)   = ()

function alg_cache(alg::QNDF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  ηold = one(uToltype)
  uprev2 = u
  uprev3 = u
  tprev2 = t
  dtₙ₋₁ = t

  D = fill(zero(typeof(t)), 1, 1)
  temp_D = fill(zero(typeof(t)), 1, 1)
  D2 = fill(zero(typeof(t)), 1, 1)
  R = fill(zero(typeof(t)), 1, 1)
  U = fill(zero(typeof(t)), 1, 1)

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

  QNDF1ConstantCache(uf,ηold,κ,tol,10000,eulercache,D,temp_D,D2,R,U,uprev2,uprev3,tprev2,dtₙ₋₁)
end

function alg_cache(alg::QNDF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u))
  W = similar(J)
  zprev = similar(u,indices(u))
  zₙ₋₁ = similar(u,indices(u)); z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zero(rate_prototype)
  fsalfirstprev = zero(rate_prototype)
  k = zero(rate_prototype)

  D = Vector{typeof(u)}(1)
  R = Vector{typeof(u)}(1)
  U = Vector{typeof(u)}(1)
  D2 = Vector{typeof(u)}(1)
  temp_D = Vector{typeof(u)}(1)

  D[1,1] = similar(u)
  R[1,1] = similar(u)
  U[1,1] = similar(u)
  D2[1,1] = similar(u)
  temp_D[1,1] = similar(u)
  temp_u = similar(u)

  tmp = similar(u); b = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))

  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  uToltype = real(uBottomEltypeNoUnits)
  utilde = similar(u,indices(u))

  uprev2 = similar(u)

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

  eulercache = ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000)

  dtₙ₋₁ = one(dt)
  QNDF1Cache(uprev2,du1,fsalfirst,fsalfirstprev,k,z,zₙ₋₁,dz,b,D,temp_D,temp_u,D2,R,U,tmp,atmp,utilde,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,eulercache,dtₙ₋₁)
end
