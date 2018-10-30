DiffEqBase.@def iipnlcachefields begin
  nlcache = alg.nlsolve.cache
  @unpack κ,tol,max_iter,min_iter,new_W = nlcache
  z = similar(u)
  dz = similar(u); tmp = similar(u); b = similar(u)
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  uToltype = real(uBottomEltypeNoUnits)
  ηold = one(uToltype)

  nf = f isa SplitFunction && alg isa SplitAlgorithms ? f.f1 : f
  islin = (f isa ODEFunction && islinear(f.f)) || (f isa SplitFunction && islinear(f.f1.f))
  # check if `nf` is linear
  if islin && alg.nlsolve isa NLNewton
    # get the operator
    J = nf.f
    W = WOperator(f.mass_matrix, dt, J)
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = alg.linsolve(Val{:init},nf,u)
    z₊ = z
  elseif alg.nlsolve isa NLNewton
    if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
      W = WOperator(f, dt)
      J = nothing # is J = W.J better?
    else
      J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
      W = similar(J)
    end
    du1 = zero(rate_prototype)
    # if the algorithm specializes on split problems the use `nf`
    uf = DiffEqDiffTools.UJacobianWrapper(nf,t,p)
    jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
    linsolve = alg.linsolve(Val{:init},uf,u)
    z₊ = z
  elseif typeof(alg.nlsolve) <: NLFunctional
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
  _nlsolve = alg.nlsolve
end
DiffEqBase.@def oopnlcachefields begin
  nlcache = alg.nlsolve.cache
  @unpack κ,tol,max_iter,min_iter,new_W = nlcache
  z = uprev
  nf = f isa SplitFunction && alg isa SplitAlgorithms ? f.f1 : f
  if alg.nlsolve isa NLNewton
    # only use `nf` if the algorithm specializes on split eqs
    uf = DiffEqDiffTools.UDerivativeWrapper(nf,t,p)
  else
    uf = nothing
  end
  islin = (f isa ODEFunction && islinear(f.f)) || (f isa SplitFunction && islinear(f.f1.f))
  if (islin || DiffEqBase.has_jac(f)) && typeof(alg.nlsolve) <: NLNewton
    # get the operator
    J = islin ? nf.f : f.jac(uprev, p, t)
    if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
      J = DiffEqArrayOperator(J)
    end
    W = WOperator(f.mass_matrix, dt, J)
  else
    W = typeof(u) <: Number ? u : Matrix{uEltypeNoUnits}(undef, 0, 0) # uEltype?
  end
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
  z₊,dz,tmp,b,k = z,z,z,z,rate_prototype
  _nlsolve = oop_nlsolver(alg.nlsolve)
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
  atmp = similar(u,uEltypeNoUnits)
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve)
end

mutable struct ImplicitEulerConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))
  ImplicitEulerConstantCache(uf,nlsolve)
end

mutable struct ImplicitMidpointConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1//2,ηold,z₊,dz,tmp,b,k))
  ImplicitMidpointConstantCache(uf,nlsolve)
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
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1//2,ηold,z₊,dz,tmp,b,k))
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
  @oopnlcachefields
  uprev3 = u
  tprev2 = t
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1,ηold,z₊,dz,tmp,b,k))
  TrapezoidConstantCache(uf,uprev3,tprev2,nlsolve)
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
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits)
  uprev3 = similar(u)
  tprev2 = t
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//2,1,ηold,z₊,dz,tmp,b,k))
  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,uprev3,tprev2,nlsolve)
end

mutable struct TRBDF2ConstantCache{F,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  tab = TRBDF2Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.d,tab.γ,ηold,z₊,dz,tmp,b,k))
  TRBDF2ConstantCache(uf,nlsolve,tab)
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
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits); zprev = similar(u);
  zᵧ = similar(u)
  tab = TRBDF2Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.d,tab.γ,ηold,z₊,dz,tmp,b,k))
  TRBDF2Cache(u,uprev,du1,fsalfirst,k,zprev,zᵧ,z,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolve,tab)
end

mutable struct SDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))

  SDIRK2ConstantCache(uf,nlsolve)
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
  @iipnlcachefields
  z₁ = similar(u); z₂ = z;
  atmp = similar(u,uEltypeNoUnits)
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1,1,ηold,z₊,dz,tmp,b,k))

  SDIRK2Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolve)
end

mutable struct SSPSDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//4,1//1,ηold,z₊,dz,tmp,b,k))

  SSPSDIRK2ConstantCache(uf,nlsolve)
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
  @iipnlcachefields
  z₁ = similar(u); z₂ = z;
  atmp = similar(u,uEltypeNoUnits)
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,1//4,1//1,ηold,z₊,dz,tmp,b,k))

  SSPSDIRK2Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,J,
                W,uf,jac_config,linsolve,nlsolve)
end

mutable struct Kvaerno3ConstantCache{UF,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::UF
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  tab = Kvaerno3Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,2tab.γ,ηold,z₊,dz,tmp,b,k))

  Kvaerno3ConstantCache(uf,nlsolve,tab)
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
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = z;
  atmp = similar(u,uEltypeNoUnits)
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,2tab.γ,ηold,z₊,dz,tmp,b,k))

  Kvaerno3Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolve,tab)
end

mutable struct Cash4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  tab = Cash4Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.γ,ηold,z₊,dz,tmp,b,k))

  Cash4ConstantCache(uf,nlsolve,tab)
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
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = similar(u)
  z₅ = z
  atmp = similar(u,uEltypeNoUnits)
  tab = Cash4Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.γ,ηold,z₊,dz,tmp,b,k))

  Cash4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolve,tab)
end

mutable struct Hairer4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  if typeof(alg) <: Hairer4
    tab = Hairer4Tableau(uToltype,real(tTypeNoUnits))
  else
    tab = Hairer42Tableau(uToltype,real(tTypeNoUnits))
  end
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.γ,ηold,z₊,dz,tmp,b,k))

  Hairer4ConstantCache(uf,nlsolve,tab)
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
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = similar(u)
  z₅ = z
  dz = similar(u)
  atmp = similar(u,uEltypeNoUnits)

  if typeof(alg) <: Hairer4
    tab = Hairer4Tableau(uToltype,real(tTypeNoUnits))
  else # Hairer42
    tab = Hairer42Tableau(uToltype,real(tTypeNoUnits))
  end
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.γ,ηold,z₊,dz,tmp,b,k))

  Hairer4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
               W,uf,jac_config,linsolve,nlsolve,tab)
end
