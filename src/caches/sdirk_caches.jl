mutable struct ImplicitEulerCache{uType,rateType,uNoUnitsType,J,W,UF,JC,uToltype,F} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f)
    J = deepcopy(f.jac_prototype)
    if isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
      W = WOperator(f.mass_matrix, dt, J)
    else
      W = similar(J)
    end
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
  z = similar(u,axes(u))
  dz = similar(u,axes(u)); tmp = similar(u,axes(u)); b = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  atmp = similar(u,uEltypeNoUnits,axes(u))

  uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
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
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end

mutable struct ImplicitEulerConstantCache{F,uToltype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
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

  ImplicitEulerConstantCache(uf,ηold,κ,tol,100000)
end

mutable struct ImplicitMidpointConstantCache{F,uToltype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
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

  ImplicitMidpointConstantCache(uf,ηold,κ,tol,10000)
end

mutable struct ImplicitMidpointCache{uType,rateType,J,UF,JC,uToltype,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
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

u_cache(c::ImplicitMidpointCache)    = (c.z,c.dz)
du_cache(c::ImplicitMidpointCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u,axes(u))
  dz = similar(u,axes(u))
  tmp = similar(u); b = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)

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

  ImplicitMidpointCache(u,uprev,du1,fsalfirst,k,z,dz,b,tmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end

mutable struct TrapezoidConstantCache{F,uToltype,uType,tType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uToltype = real(uBottomEltypeNoUnits)
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  ηold = one(uToltype)
  uprev3 = u
  tprev2 = t

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

  TrapezoidConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct TrapezoidCache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,tType,F} <: OrdinaryDiffEqMutableCache
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
  uprev3::uType
  tprev2::tType
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u,axes(u))
  dz = similar(u,axes(u))
  tmp = similar(u); b = similar(u,axes(u));
  atmp = similar(u,uEltypeNoUnits,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)

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

  uprev3 = similar(u)
  tprev2 = t

  ηold = one(uToltype)

  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct TRBDF2ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  tab = TRBDF2Tableau(uToltype,real(tTypeNoUnits))

  TRBDF2ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct TRBDF2Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
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
  tab::Tab
end

u_cache(c::TRBDF2Cache)    = (c.zprev,c.zᵧ,c.z,c.dz)
du_cache(c::TRBDF2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  zprev = similar(u,axes(u));
  zᵧ = similar(u,axes(u)); z = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
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

  tab = TRBDF2Tableau(uToltype,real(tTypeNoUnits))

  ηold = one(uToltype)

  TRBDF2Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,zprev,zᵧ,z,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct SDIRK2ConstantCache{F,uToltype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  SDIRK2ConstantCache(uf,ηold,κ,tol,10000)
end

mutable struct SDIRK2Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
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

u_cache(c::SDIRK2Cache)    = (c.z₁,c.z₂,c.dz)
du_cache(c::SDIRK2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u));
  z₂ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
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

  SDIRK2Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end

mutable struct SSPSDIRK2ConstantCache{F,uToltype} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  uToltype = real(uBottomEltypeNoUnits)
  ηold = one(uToltype)
  uprev3 = u
  tprev2 = t

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

  SSPSDIRK2ConstantCache(uf,ηold,κ,tol,10000)
end

mutable struct SSPSDIRK2Cache{uType,rateType,J,UF,JC,uToltype,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  dz::uType
  b::uType
  tmp::uType
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

u_cache(c::SSPSDIRK2Cache)    = (c.z₁,c.z₂,c.dz)
du_cache(c::SSPSDIRK2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u));
  z₂ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
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

  SSPSDIRK2Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end

mutable struct Kvaerno3ConstantCache{UF,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::UF
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  tab = Kvaerno3Tableau(uToltype,real(tTypeNoUnits))

  Kvaerno3ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno3Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
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
  tab::Tab
end

u_cache(c::Kvaerno3Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.dz)
du_cache(c::Kvaerno3Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u));
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
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

  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uToltype)

  Kvaerno3Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct Cash4ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  tab = Cash4Tableau(uToltype,real(tTypeNoUnits))

  Cash4ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct Cash4Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
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
  tab::Tab
end

u_cache(c::Cash4Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.dz)
du_cache(c::Cash4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u));
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  z₅ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
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

  tab = Cash4Tableau(uToltype,real(tTypeNoUnits))

  ηold = one(uToltype)

  Cash4Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct Hairer4ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  if typeof(alg) <: Hairer4
    tab = Hairer4Tableau(uToltype,real(tTypeNoUnits))
  else
    tab = Hairer42Tableau(uToltype,real(tTypeNoUnits))
  end

  Hairer4ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct Hairer4Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
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
  tab::Tab
end

u_cache(c::Hairer4Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.dz)
du_cache(c::Hairer4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u));
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  z₅ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
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

  if typeof(alg) <: Hairer4
    tab = Hairer4Tableau(uToltype,real(tTypeNoUnits))
  else # Hairer42
    tab = Hairer42Tableau(uToltype,real(tTypeNoUnits))
  end

  ηold = one(uToltype)

  Hairer4Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end
