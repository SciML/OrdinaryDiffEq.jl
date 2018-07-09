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
  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  zprev = similar(u,indices(u))
  zₙ₋₁ = similar(u,indices(u)); z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zeros(rate_prototype)
  fsalfirstprev = zeros(rate_prototype)
  k = zeros(rate_prototype)
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

  eulercache = ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000)

  dtₙ₋₁ = one(dt)
  ABDF2Cache(u,uprev,uprev2,du1,fsalfirst,fsalfirstprev,k,z,zₙ₋₁,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,eulercache,dtₙ₋₁)
end

# QNDF1

mutable struct QNDF1ConstantCache{F,uToltype,coefType,coefType1,dtType,uType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  D::coefType
  D2::coefType1
  R::coefType
  U::coefType
  uprev2::uType
  dtₙ₋₁::dtType
end

mutable struct QNDF1Cache{uType,rateType,coefType,uNoUnitsType,J,UF,JC,uToltype,F,dtType} <: OrdinaryDiffEqMutableCache
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  b::uType
  D::coefType
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
  dtₙ₋₁::dtType
end

u_cache(c::QNDF1Cache)    = ()
du_cache(c::QNDF1Cache)   = ()

function alg_cache(alg::QNDF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  ηold = one(uToltype)
  uprev2 = u
  dtₙ₋₁ = t

  D = zeros(typeof(t), 1, 1)
  D2 = zeros(typeof(t), 1, 2)
  R = zeros(typeof(t), 1, 1)
  U = zeros(typeof(t), 1, 1)

  U!(1,U,false)

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

  QNDF1ConstantCache(uf,ηold,κ,tol,10000,D,D2,R,U,uprev2,dtₙ₋₁)
end

function alg_cache(alg::QNDF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u))
  W = similar(J)
  z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)

  D = Vector{typeof(u)}(1)
  R = Vector{typeof(u)}(1)
  U = Vector{typeof(u)}(1)
  D2 = Vector{typeof(u)}(2)

  D[1] = similar(u)
  R[1] = similar(u)
  U[1] = similar(u)
  D2[1] = similar(u)
  D2[2] = similar(u)

  U!(1,U,true)

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
  dtₙ₋₁ = one(dt)

  QNDF1Cache(uprev2,du1,fsalfirst,k,z,dz,b,D,D2,R,U,tmp,atmp,utilde,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,dtₙ₋₁)
end

# QNDF2

mutable struct QNDF2ConstantCache{F,uToltype,coefType,coefType1,uType,dtType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  eulercache::ImplicitEulerConstantCache
  D::coefType
  D2::coefType
  R::coefType1
  U::coefType1
  uprev2::uType
  uprev3::uType
  dtₙ₋₁::dtType
  dtₙ₋₂::dtType
  pass1::Bool
  pass2::Bool
  pdt::dtType
end

function alg_cache(alg::QNDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  ηold = one(uToltype)
  uprev2 = u
  uprev3 = u
  dtₙ₋₁ = zero(t)
  dtₙ₋₂ = one(t)
  pdt = t

  D = zeros(typeof(u),1,2)
  D2 = zeros(typeof(u),1,3)
  R = zeros(typeof(u),2,2)
  U = zeros(typeof(u),2,2)

  U!(2,U,false)

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

  QNDF2ConstantCache(uf,ηold,κ,tol,10000,eulercache,D,D2,R,U,uprev2,uprev3,dtₙ₋₁,dtₙ₋₂,true,true,pdt)
end
