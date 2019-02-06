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
  γ, c = 1, 1
  @iipnlsolve

  atmp = similar(u,uEltypeNoUnits)

  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolve)
end

mutable struct ImplicitEulerConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1, 1
  @oopnlsolve
  ImplicitEulerConstantCache(uf,nlsolve)
end

mutable struct ImplicitMidpointConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
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
  γ, c = 1//2, 1
  @oopnlsolve

  uprev3 = u
  tprev2 = t

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
  γ, c = 1//2, 1
  @iipnlsolve

  uprev3 = similar(u)
  tprev2 = t
  atmp = similar(u,uEltypeNoUnits)

  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,uprev3,tprev2,nlsolve)
end

mutable struct TRBDF2ConstantCache{F,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = TRBDF2Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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
  tab = TRBDF2Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  @iipnlsolve

  atmp = similar(u,uEltypeNoUnits); zprev = similar(u); zᵧ = similar(u)

  TRBDF2Cache(u,uprev,du1,fsalfirst,k,zprev,zᵧ,z,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolve,tab)
end

mutable struct SDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
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
  γ, c = 1, 1
  @iipnlsolve

  z₁ = similar(u); z₂ = z
  atmp = similar(u,uEltypeNoUnits)

  SDIRK2Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolve)
end

mutable struct SSPSDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
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
  γ, c = 1//4, 1//1
  @iipnlsolve

  z₁ = similar(u); z₂ = z
  atmp = similar(u,uEltypeNoUnits)

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
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  @iipnlsolve

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = z
  atmp = similar(u,uEltypeNoUnits)

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
  tab = Cash4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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
  tab = Cash4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  @iipnlsolve

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u); z₅ = z
  atmp = similar(u,uEltypeNoUnits)

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
  if alg isa Hairer4
    tab = Hairer4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  else
    tab = Hairer42Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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
  if alg isa Hairer4
    tab = Hairer4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  else # Hairer42
    tab = Hairer42Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  end
  γ, c = tab.γ, tab.γ
  @iipnlsolve

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u); z₅ = z
  atmp = similar(u,uEltypeNoUnits)

  Hairer4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
               W,uf,jac_config,linsolve,nlsolve,tab)
end
