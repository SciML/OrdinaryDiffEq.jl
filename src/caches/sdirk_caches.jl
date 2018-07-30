DiffEqBase.@def iipnlcachefields begin
  nlcache = alg.nonlinsolve.cache
  @unpack κ,tol,max_iter,min_iter,new_W = nlcache
  z = similar(u,axes(u))
  dz = similar(u,axes(u)); tmp = similar(u,axes(u)); b = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  uToltype = real(uBottomEltypeNoUnits)
  ηold = one(uToltype)

  if typeof(alg.nonlinsolve) <: NLNewton
    if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
      W = WOperator(f, dt)
      J = nothing # is J = W.J better?
    else
      J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
      W = similar(J)
    end
    du1 = zero(rate_prototype)
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
    jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
    linsolve = alg.linsolve(Val{:init},uf,u)
    z₊ = z
  elseif typeof(alg.nonlinsolve) <: NLFunctional
    J = nothing
    W = nothing
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = nothing
    z₊ = similar(z)
  end

  if κ != nothing
    κ = uToltype(nlcache.κ)
  else
    κ = uToltype(1//100)
  end
  if tol != nothing
    tol = uToltype(nlcache.tol)
  else
    tol = uToltype(min(0.03,first(reltol)^(0.5)))
  end
end
DiffEqBase.@def oopnlcachefields begin
  nlcache = alg.nonlinsolve.cache
  @unpack κ,tol,max_iter,min_iter,new_W = nlcache
  z = uprev
  if typeof(alg.nonlinsolve) <: NLNewton
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  else
    uf = nothing
  end
  W = typeof(u) <: Number ? u : Matrix{uEltypeNoUnits}(undef, 0, 0) # uEltype?
  uToltype = real(uBottomEltypeNoUnits)
  ηold = one(uToltype)

  if κ != nothing
    κ = uToltype(nlcache.κ)
  else
    κ = uToltype(1//100)
  end
  if tol != nothing
    tol = uToltype(nlcache.tol)
  else
    tol = uToltype(min(0.03,first(reltol)^(0.5)))
  end
end

mutable struct ImplicitEulerCache{uType,rateType,uNoUnitsType,J,W,UF,JC,F,N} <: OrdinaryDiffEqMutableCache
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
  nlsolve::N
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits,axes(u))
  nlsolve = typeof(alg.nonlinsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve)
end

mutable struct ImplicitEulerConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(alg.nonlinsolve)(NLSolverCache(κ,tol,min_iter,max_iter,100000,new_W,z,W,1,1,ηold,z,z,z,z,rate_prototype))
  ImplicitEulerConstantCache(uf,oop_nlsolver(nlsolve))
end

mutable struct ImplicitMidpointConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(alg.nonlinsolve)(NLSolverCache(κ,tol,min_iter,max_iter,100000,new_W,z,W,1//2,1//2,ηold,z,z,z,z,rate_prototype))
  ImplicitMidpointConstantCache(uf,oop_nlsolver(nlsolve))
end

mutable struct ImplicitMidpointCache{uType,rateType,J,W,UF,JC,F,N} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

u_cache(c::ImplicitMidpointCache)    = (c.z,c.dz)
du_cache(c::ImplicitMidpointCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  nlsolve = typeof(alg.nonlinsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1//2,ηold,z₊,dz,tmp,b,k))
  ImplicitMidpointCache(u,uprev,du1,fsalfirst,k,z,dz,b,tmp,J,W,uf,jac_config,linsolve,nlsolve)
end

mutable struct TrapezoidConstantCache{F,uType,tType,N} <: OrdinaryDiffEqConstantCache
  uf::F
  uprev3::uType
  tprev2::tType
  nlsolve::N
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

mutable struct TrapezoidCache{uType,rateType,uNoUnitsType,J,W,UF,JC,tType,F,N} <: OrdinaryDiffEqMutableCache
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
  uprev3::uType
  tprev2::tType
  nlsolve::N
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
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

mutable struct TRBDF2ConstantCache{F,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
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

mutable struct TRBDF2Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,Tab,F,N} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

u_cache(c::TRBDF2Cache)    = (c.zprev,c.zᵧ,c.z,c.dz)
du_cache(c::TRBDF2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
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

  TRBDF2Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,zprev,zᵧ,z,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct SDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
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

mutable struct SDIRK2Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,F,N} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

u_cache(c::SDIRK2Cache)    = (c.z₁,c.z₂,c.dz)
du_cache(c::SDIRK2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
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

  SDIRK2Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uToltype,typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end

mutable struct SSPSDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
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

mutable struct SSPSDIRK2Cache{uType,rateType,J,W,UF,JC,F,N} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

u_cache(c::SSPSDIRK2Cache)    = (c.z₁,c.z₂,c.dz)
du_cache(c::SSPSDIRK2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
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

  SSPSDIRK2Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uToltype,typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end

mutable struct Kvaerno3ConstantCache{UF,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::UF
  nlsolve::N
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

mutable struct Kvaerno3Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,Tab,F,N} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

u_cache(c::Kvaerno3Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.dz)
du_cache(c::Kvaerno3Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
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

  Kvaerno3Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct Cash4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
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

mutable struct Cash4Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,N,Tab,F} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

u_cache(c::Cash4Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.dz)
du_cache(c::Cash4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
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

  Cash4Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct Hairer4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
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

mutable struct Hairer4Cache{uType,rateType,uNoUnitsType,J,W,UF,JC,Tab,F,N} <: OrdinaryDiffEqMutableCache
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
  W::W
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

u_cache(c::Hairer4Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.dz)
du_cache(c::Hairer4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, dt)
    J = nothing # is J = W.J better?
  else
    J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
    W = similar(J)
  end
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

  Hairer4Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end
