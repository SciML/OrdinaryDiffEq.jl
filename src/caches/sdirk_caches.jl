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
    W = WOperator(f.mass_matrix, dt, J, true)
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = alg.linsolve(Val{:init},nf,u)
    z₊ = z
  elseif alg.nlsolve isa NLNewton
    if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype != nothing
      W = WOperator(f, dt, inplace)
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
    W = WOperator(f.mass_matrix, dt, J, false)
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
end

DiffEqBase.@def iipnlsolve begin
  _nlsolve = typeof(alg.nlsolve).name.wrapper
  nlcache = NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,γ,c,ηold,z₊,dz,tmp,b,k)
  nlsolve = _nlsolve{true, typeof(nlcache)}(nlcache)
end
DiffEqBase.@def oopnlsolve begin
  _nlsolve = typeof(alg.nlsolve).name.wrapper
  nlcache = NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,γ,c,ηold,z₊,dz,tmp,b,k)
  nlsolve = _nlsolve{false, typeof(nlcache)}(nlcache)
end

abstract type SDIRKMutableCache <: OrdinaryDiffEqMutableCache end

@cache mutable struct ImplicitEulerCache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits)
  γ, c = 1, 1
  @iipnlsolve
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve)
end

mutable struct ImplicitEulerConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  γ, c = 1, 1
  @oopnlsolve
  ImplicitEulerConstantCache(uf,nlsolve)
end

mutable struct ImplicitMidpointConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  γ, c = 1//2, 1//2
  @oopnlsolve
  ImplicitMidpointConstantCache(uf,nlsolve)
end

@cache mutable struct ImplicitMidpointCache{uType,rateType,JType,WType,UF,JC,F,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  b::uType
  tmp::uType
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  γ, c = 1//2, 1//2
  @iipnlsolve
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
  γ, c = 1//2, 1
  @oopnlsolve
  TrapezoidConstantCache(uf,uprev3,tprev2,nlsolve)
end

@cache mutable struct TrapezoidCache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,tType,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  uprev3::uType
  tprev2::tType
  nlsolve::N
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits)
  uprev3 = similar(u)
  tprev2 = t
  γ, c = 1//2, 1
  @iipnlsolve
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
  γ, c = tab.d, tab.γ
  @oopnlsolve
  TRBDF2ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct TRBDF2Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,Tab,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  atmp = similar(u,uEltypeNoUnits); zprev = similar(u);
  zᵧ = similar(u)
  tab = TRBDF2Tableau(uToltype,real(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  @iipnlsolve
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
  γ, c = 1, 1
  @oopnlsolve
  SDIRK2ConstantCache(uf,nlsolve)
end

@cache mutable struct SDIRK2Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = z;
  atmp = similar(u,uEltypeNoUnits)
  γ, c = 1, 1
  @iipnlsolve
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
  γ, c = 1//4, 1//1
  @oopnlsolve
  SSPSDIRK2ConstantCache(uf,nlsolve)
end

@cache mutable struct SSPSDIRK2Cache{uType,rateType,JType,WType,UF,JC,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = z;
  atmp = similar(u,uEltypeNoUnits)
  γ, c = 1//4, 1//1
  @iipnlsolve
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
  γ, c = tab.γ, 2tab.γ
  @oopnlsolve
  Kvaerno3ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct Kvaerno3Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,Tab,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = z;
  atmp = similar(u,uEltypeNoUnits)
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  @iipnlsolve

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
  γ, c = tab.γ,tab.γ
  @oopnlsolve
  Cash4ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct Cash4Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = similar(u)
  z₅ = z
  atmp = similar(u,uEltypeNoUnits)
  tab = Cash4Tableau(uToltype,real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  @iipnlsolve
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
  γ, c = tab.γ, tab.γ
  @oopnlsolve

  Hairer4ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct Hairer4Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,Tab,F,N} <: SDIRKMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolve::N
  tab::Tab
end

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
  γ, c = tab.γ, tab.γ
  @iipnlsolve

  Hairer4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
               W,uf,jac_config,linsolve,nlsolve,tab)
end
