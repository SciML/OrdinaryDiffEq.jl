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

mutable struct QNDF1Cache{uType,rateType,coefType,coefType1,coefType2,uNoUnitsType,J,UF,JC,uToltype,F,dtType} <: OrdinaryDiffEqMutableCache
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  b::uType
  D::coefType1
  D2::coefType2
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

  D = fill(zero(typeof(u)), 1, 1)
  D2 = fill(zero(typeof(u)), 1, 2)
  R = fill(zero(typeof(t)), 1, 1)
  U = fill(zero(typeof(t)), 1, 1)

  U!(1,U)

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
  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u))
  W = similar(J)
  z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)

  D = Array{typeof(u)}(undef, 1, 1)
  D2 = Array{typeof(u)}(undef, 1, 2)
  R = fill(zero(typeof(t)), 1, 1)
  U = fill(zero(typeof(t)), 1, 1)
 
  D[1] = similar(u)
  D2[1] = similar(u); D2[2] = similar(u)

  U!(1,U)

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
  D::coefType
  D2::coefType
  R::coefType1
  U::coefType1
  uprev2::uType
  uprev3::uType
  dtₙ₋₁::dtType
  dtₙ₋₂::dtType
end

mutable struct QNDF2Cache{uType,rateType,coefType,coefType1,coefType2,uNoUnitsType,J,UF,JC,uToltype,F,dtType} <: OrdinaryDiffEqMutableCache
  uprev2::uType
  uprev3::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  b::uType
  D::coefType1
  D2::coefType2
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
  dtₙ₋₂::dtType
end

u_cache(c::QNDF2Cache)  = ()
du_cache(c::QNDF2Cache) = ()

function alg_cache(alg::QNDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  ηold = one(uToltype)
  uprev2 = u
  uprev3 = u
  dtₙ₋₁ = zero(t)
  dtₙ₋₂ = zero(t)

  D = fill(zero(typeof(u)), 1, 2)
  D2 = fill(zero(typeof(u)), 1, 3)
  R = fill(zero(typeof(t)), 2, 2)
  U = fill(zero(typeof(t)), 2, 2)

  U!(2,U)

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

  QNDF2ConstantCache(uf,ηold,κ,tol,10000,D,D2,R,U,uprev2,uprev3,dtₙ₋₁,dtₙ₋₂)
end

function alg_cache(alg::QNDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u))
  W = similar(J)
  z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)

  D = Array{typeof(u)}(undef, 1, 2)
  D2 = Array{typeof(u)}(undef, 1, 3)
  R = fill(zero(typeof(t)), 2, 2)
  U = fill(zero(typeof(t)), 2, 2)
  
  D[1] = similar(u); D[2] = similar(u)
  D2[1] = similar(u);  D2[2] = similar(u); D2[3] = similar(u)

  U!(2,U)

  tmp = similar(u); b = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))

  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  uToltype = real(uBottomEltypeNoUnits)
  utilde = similar(u,indices(u))

  uprev2 = similar(u)
  uprev3 = similar(u)

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
  dtₙ₋₁ = zero(dt)
  dtₙ₋₂ = zero(dt)

  QNDF2Cache(uprev2,uprev3,du1,fsalfirst,k,z,dz,b,D,D2,R,U,tmp,atmp,utilde,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,dtₙ₋₁,dtₙ₋₂)
end

mutable struct QNDFConstantCache{F,uToltype,coefType,coefType1,uType,dtType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  D::coefType
  D2::coefType2
  R::coefType1
  U::coefType1
  k::Int64
  max_order::Int64
  udiff::uType
  dts::dtType
end


function alg_cache(alg::QNDF,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  ηold = one(uToltype)
  udiff = fill(zero(typeof(u)), 1, 6)
  dts = fill(zero(typeof(u)), 1, 6)

  D = fill(zero(typeof(u)), 1, 6)
  D2 = fill(zero(typeof(u)), 6, 6)
  R = fill(zero(typeof(t)), 5, 5)
  U = fill(zero(typeof(t)), 5, 5)

  max_order = 5

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

  QNDFConstantCache(uf,ηold,κ,tol,10000,D,D2,R,U,1,max_order,udiff,dts)
end
