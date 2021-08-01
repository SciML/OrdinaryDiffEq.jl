abstract type SDIRKMutableCache <: OrdinaryDiffEqMutableCache end

@cache mutable struct ImplicitEulerCache{uType,rateType,uNoUnitsType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  atmp::uNoUnitsType
  nlsolver::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  atmp = similar(u,uEltypeNoUnits)

  ImplicitEulerCache(u,uprev,uprev2,fsalfirst,atmp,nlsolver)
end

mutable struct ImplicitEulerConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  ImplicitEulerConstantCache(nlsolver)
end

mutable struct ImplicitMidpointConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//2, 1//2
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  ImplicitMidpointConstantCache(nlsolver)
end

@cache mutable struct ImplicitMidpointCache{uType,rateType,N} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  nlsolver::N
end

function alg_cache(alg::ImplicitMidpoint,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//2, 1//2
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)
  ImplicitMidpointCache(u,uprev,fsalfirst,nlsolver)
end

mutable struct TrapezoidConstantCache{uType,tType,N} <: OrdinaryDiffEqConstantCache
  uprev3::uType
  tprev2::tType
  nlsolver::N
end

function alg_cache(alg::Trapezoid,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//2, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

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

function alg_cache(alg::Trapezoid,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//2, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
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

function alg_cache(alg::TRBDF2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = TRBDF2Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
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

function alg_cache(alg::TRBDF2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = TRBDF2Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.d, tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  atmp = similar(u,uEltypeNoUnits); zprev = zero(u); zᵧ = zero(u)

  TRBDF2Cache(u,uprev,fsalfirst,zprev,zᵧ,atmp,nlsolver,tab)
end

mutable struct SDIRK2ConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::SDIRK2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
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

function alg_cache(alg::SDIRK2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SDIRK2Cache(u,uprev,fsalfirst,z₁,z₂,atmp,nlsolver)
end

struct SDIRK22ConstantCache{uType,tType,N,Tab} <: OrdinaryDiffEqConstantCache
  uprev3::uType
  tprev2::tType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SDIRK22,u,rate_prototype,::Type{uEltypeNoUnits},::Type{tTypeNoUnits},::Type{uBottomEltypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SDIRK22Tableau(constvalue(uBottomEltypeNoUnits))
  uprev3 = u
  tprev2 = t
  γ, c = 1, 1

  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  SDIRK22ConstantCache(uprev3,tprev2,nlsolver)
end

@cache mutable struct SDIRK22Cache{uType,rateType,uNoUnitsType,tType,N,Tab} <: SDIRKMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  atmp::uNoUnitsType
  uprev3::uType
  tprev2::tType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SDIRK22,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SDIRK22Tableau(constvalue(uBottomEltypeNoUnits))
  γ, c = 1, 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  uprev3 = zero(u)
  tprev2 = t
  atmp = similar(u,uEltypeNoUnits)

  SDIRK22(u,uprev,uprev2,fsalfirst,atmp,uprev3,tprev2,nlsolver,tab)
end

mutable struct SSPSDIRK2ConstantCache{N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//4, 1//1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
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

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = 1//4, 1//1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SSPSDIRK2Cache(u,uprev,fsalfirst,z₁,z₂,nlsolver)
end

mutable struct Kvaerno3ConstantCache{Tab,N} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = Kvaerno3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
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

function alg_cache(alg::Kvaerno3,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = Kvaerno3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ, 2tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Kvaerno3Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,atmp,nlsolver,tab)
end

mutable struct Cash4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = Cash4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
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

function alg_cache(alg::Cash4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = Cash4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u); z₅ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Cash4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,atmp,nlsolver,tab)
end

mutable struct SFSDIRK4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  SFSDIRK4ConstantCache(nlsolver,tab)
end

@cache mutable struct SFSDIRK4Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
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

function alg_cache(alg::SFSDIRK4,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u); z₅ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SFSDIRK4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,atmp,nlsolver,tab)
end

mutable struct SFSDIRK5ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  SFSDIRK5ConstantCache(nlsolver,tab)
end

@cache mutable struct SFSDIRK5Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  z₆::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u);z₅ = zero(u); z₆ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SFSDIRK5Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,atmp,nlsolver,tab)
end

mutable struct SFSDIRK6ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK6,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK6Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  SFSDIRK6ConstantCache(nlsolver,tab)
end

@cache mutable struct SFSDIRK6Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  z₆::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK6,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK6Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u);z₅ = zero(u); z₆ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SFSDIRK6Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,atmp,nlsolver,tab)
end

mutable struct SFSDIRK7ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK7,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK7Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  SFSDIRK7ConstantCache(nlsolver,tab)
end

@cache mutable struct SFSDIRK7Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  z₆::uType
  z₇::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK7,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK7Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u);z₅ = zero(u);z₆ = zero(u); z₇ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SFSDIRK7Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,z₇,atmp,nlsolver,tab)
end

mutable struct SFSDIRK8ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK8,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK8Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  SFSDIRK8ConstantCache(nlsolver,tab)
end

@cache mutable struct SFSDIRK8Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  z₆::uType
  z₇::uType
  z₈::uType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SFSDIRK8,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = SFSDIRK8Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u);z₅ = zero(u);z₆ = zero(u);z₇ = zero(u); z₈ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  SFSDIRK8Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,atmp,nlsolver,tab)
end

mutable struct Hairer4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  if alg isa Hairer4
    tab = Hairer4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  else
    tab = Hairer42Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  end
  γ, c = tab.γ, tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
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

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  if alg isa Hairer4
    tab = Hairer4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  else # Hairer42
    tab = Hairer42Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  end
  γ, c = tab.γ, tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u); z₅ = nlsolver.z
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

function alg_cache(alg::ESDIRK54I8L2SA,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = ESDIRK54I8L2SATableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ, tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
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

function alg_cache(alg::ESDIRK54I8L2SA,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  tab = ESDIRK54I8L2SATableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  γ, c = tab.γ,tab.γ
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  ESDIRK54I8L2SAConstantCache(nlsolver,tab)
end
