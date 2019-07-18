abstract type SDIRKMutableCache <: OrdinaryDiffEqMutableCache end

@cache mutable struct ImplicitEulerCache{uType,rateType,uNoUnitsType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  atmp::uNoUnitsType
  nlsolver::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  atmp = similar(u,uEltypeNoUnits)

  ImplicitEulerCache(u,uprev,uprev2,fsalfirst,atmp,nlsolver)
end

mutable struct ImplicitEulerConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1, 1
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  ImplicitEulerConstantCache(nlsolver)
end

mutable struct ImplicitMidpointConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//2, 1//2
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  ImplicitMidpointConstantCache(nlsolver)
end

@cache mutable struct ImplicitMidpointCache{uType,rateType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  nlsolver::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//2, 1//2
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)
  ImplicitMidpointCache(u,uprev,fsalfirst,nlsolver)
end

mutable struct TrapezoidConstantCache{uType,tType,N} <: OrdinaryDiffEqConstantCache
  uprev3::uType
  tprev2::tType
  nlsolver::N
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//2, 1
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)

  uprev3 = u
  tprev2 = t

  TrapezoidConstantCache(uprev3,tprev2,nlsolver)
end

@cache mutable struct TrapezoidCache{uType,rateType,uNoUnitsType,tType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  atmp::uNoUnitsType
  uprev3::uType
  tprev2::tType
  nlsolver::N
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//2, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  uprev3 = zero(u)
  tprev2 = t
  atmp = similar(u,uEltypeNoUnits)

  TrapezoidCache(u,uprev,uprev2,fsalfirst,atmp,uprev3,tprev2,nlsolver)
end

mutable struct TRBDF2ConstantCache{Tab,N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = TRBDF2Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  TRBDF2ConstantCache(nlsolver,tab)
end

@cache mutable struct TRBDF2Cache{uType,rateType,uNoUnitsType,Tab,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  zprev::uType
  zᵧ::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = TRBDF2Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  atmp = similar(u,uEltypeNoUnits); zprev = similar(u); zᵧ = similar(u)

  TRBDF2Cache(u,uprev,fsalfirst,zprev,zᵧ,atmp,nlsolver,tab)
end

mutable struct SDIRK2ConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1, 1
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  SDIRK2ConstantCache(nlsolver)
end

@cache mutable struct SDIRK2Cache{uType,rateType,uNoUnitsType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  atmp::uNoUnitsType
  nlsolver::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1, 1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  z₁ = similar(u); z₂ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SDIRK2Cache(u,uprev,fsalfirst,z₁,z₂,atmp,nlsolver)
end

mutable struct SSPSDIRK2ConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  γ, c = 1//4, 1//1
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  SSPSDIRK2ConstantCache(nlsolver)
end

@cache mutable struct SSPSDIRK2Cache{uType,rateType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  nlsolver::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  γ, c = 1//4, 1//1
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  z₁ = similar(u); z₂ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SSPSDIRK2Cache(u,uprev,fsalfirst,z₁,z₂,nlsolver)
end

mutable struct Kvaerno3ConstantCache{Tab,N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  Kvaerno3ConstantCache(nlsolver,tab)
end

@cache mutable struct Kvaerno3Cache{uType,rateType,uNoUnitsType,Tab,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = Kvaerno3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Kvaerno3Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,atmp,nlsolver,tab)
end

mutable struct Cash4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = Cash4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  Cash4ConstantCache(nlsolver,tab)
end

@cache mutable struct Cash4Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = Cash4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u); z₅ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Cash4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,atmp,nlsolver,tab)
end

mutable struct Hairer4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
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
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  Hairer4ConstantCache(nlsolver,tab)
end

@cache mutable struct Hairer4Cache{uType,rateType,uNoUnitsType,Tab,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  atmp::uNoUnitsType
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
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u); z₅ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Hairer4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,atmp,nlsolver,tab)
end

@cache mutable struct ESDIRK54I8L2SACache{uType,rateType,uNoUnitsType,Tab,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType; z₂::uType; z₃::uType; z₄::uType; z₅::uType; z₆::uType; z₇::uType; z₈::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::ESDIRK54I8L2SA,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = ESDIRK54I8L2SATableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.γ
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u)
  z₅ = zero(u); z₆ = zero(u); z₇ = zero(u); z₈ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  ESDIRK54I8L2SACache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,atmp,nlsolver,tab)
end

mutable struct ESDIRK54I8L2SAConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::ESDIRK54I8L2SA,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  tab = ESDIRK54I8L2SATableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  ESDIRK54I8L2SAConstantCache(nlsolver,tab)
end
