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
  nlsolver::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  atmp = similar(u,uEltypeNoUnits)

  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,nlsolver)
end

mutable struct ImplicitEulerConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1, 1
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  ImplicitEulerConstantCache(uf,nlsolver)
end

mutable struct ImplicitMidpointConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//2, 1//2
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  ImplicitMidpointConstantCache(uf,nlsolver)
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
  nlsolver::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//2, 1//2
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields
  ImplicitMidpointCache(u,uprev,du1,fsalfirst,k,z,dz,b,tmp,J,W,uf,jac_config,linsolve,nlsolver)
end

mutable struct TrapezoidConstantCache{F,uType,tType,N} <: OrdinaryDiffEqConstantCache
  uf::F
  uprev3::uType
  tprev2::tType
  nlsolver::N
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//2, 1
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields

  uprev3 = u
  tprev2 = t

  TrapezoidConstantCache(uf,uprev3,tprev2,nlsolver)
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
  nlsolver::N
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//2, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  uprev3 = zero(u)
  tprev2 = t
  atmp = similar(u,uEltypeNoUnits)

  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,b,tmp,atmp,J,W,uf,jac_config,linsolve,uprev3,tprev2,nlsolver)
end

mutable struct TRBDF2ConstantCache{F,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = TRBDF2Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  TRBDF2ConstantCache(uf,nlsolver,tab)
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
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = TRBDF2Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  atmp = similar(u,uEltypeNoUnits); zprev = similar(u); zᵧ = similar(u)

  TRBDF2Cache(u,uprev,du1,fsalfirst,k,zprev,zᵧ,z,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolver,tab)
end

mutable struct SDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1, 1
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  SDIRK2ConstantCache(uf,nlsolver)
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
  nlsolver::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  z₁ = similar(u); z₂ = z
  atmp = similar(u,uEltypeNoUnits)

  SDIRK2Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolver)
end

mutable struct SSPSDIRK2ConstantCache{F,N} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//4, 1//1
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  SSPSDIRK2ConstantCache(uf,nlsolver)
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
  nlsolver::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//4, 1//1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  z₁ = similar(u); z₂ = z
  atmp = similar(u,uEltypeNoUnits)

  SSPSDIRK2Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,dz,b,tmp,J,
                W,uf,jac_config,linsolve,nlsolver)
end

mutable struct Kvaerno3ConstantCache{UF,Tab,N} <: OrdinaryDiffEqConstantCache
  uf::UF
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  Kvaerno3ConstantCache(uf,nlsolver,tab)
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
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = z
  atmp = similar(u,uEltypeNoUnits)

  Kvaerno3Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolver,tab)
end

mutable struct Cash4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = Cash4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  Cash4ConstantCache(uf,nlsolver,tab)
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
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = Cash4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u); z₅ = z
  atmp = similar(u,uEltypeNoUnits)

  Cash4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolver,tab)
end

mutable struct Hairer4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
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
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  Hairer4ConstantCache(uf,nlsolver,tab)
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
  nlsolver::N
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
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u); z₅ = z
  atmp = similar(u,uEltypeNoUnits)

  Hairer4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
               W,uf,jac_config,linsolve,nlsolver,tab)
end

@cache mutable struct ESDIRK54I8L2SACache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,Tab,F,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType; z₂::uType; z₃::uType; z₄::uType; z₅::uType; z₆::uType; z₇::uType; z₈::uType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::ESDIRK54I8L2SA,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = ESDIRK54I8L2SATableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u)
  z₅ = zero(u); z₆ = zero(u); z₇ = zero(u); z₈ = z
  atmp = similar(u,uEltypeNoUnits)

  ESDIRK54I8L2SACache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,dz,b,tmp,atmp,J,
               W,uf,jac_config,linsolve,nlsolver,tab)
end

mutable struct ESDIRK54I8L2SAConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::ESDIRK54I8L2SA,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = ESDIRK54I8L2SATableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getoopnlsolvefields
  ESDIRK54I8L2SAConstantCache(uf,nlsolver,tab)
end
