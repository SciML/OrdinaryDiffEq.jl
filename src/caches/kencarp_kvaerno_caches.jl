@cache mutable struct KenCarp3ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = KenCarp3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  KenCarp3ConstantCache(nlsolver,tab)
end

@cache mutable struct KenCarp3Cache{uType,rateType,uNoUnitsType,N,Tab,kType} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = KenCarp3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  if typeof(f) <: SplitFunction
    k1 = similar(u); k2 = similar(u)
    k3 = similar(u); k4 = similar(u)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    uf = UJacobianWrapper(f,t,p)
  end

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  KenCarp3Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,k1,k2,k3,k4,atmp,nlsolver,tab)
end

@cache mutable struct CFNLIRK3ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::CFNLIRK3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = CFNLIRK3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  CFNLIRK3ConstantCache(nlsolver,tab)
end

@cache mutable struct CFNLIRK3Cache{uType,rateType,uNoUnitsType,N,Tab,kType} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::CFNLIRK3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = CFNLIRK3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  k1 = similar(u); k2 = similar(u)
  k3 = similar(u); k4 = similar(u)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  CFNLIRK3Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,k1,k2,k3,k4,atmp,nlsolver,tab)
end

@cache mutable struct CFNLIRK4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::CFNLIRK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = CFNLIRK4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  CFNLIRK4ConstantCache(nlsolver,tab)
end

@cache mutable struct CFNLIRK4Cache{uType,rateType,uNoUnitsType,N,Tab,kType} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  z₆::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  k5::kType
  k6::kType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::CFNLIRK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = CFNLIRK4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  k1 = similar(u); k2 = similar(u)
  k3 = similar(u); k4 = similar(u)
  k5 = similar(u); k6 = similar(u)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u);z₄ = zero(u); z₅ = zero(u); z₆ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  CFNLIRK4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,k1,k2,k3,k4,k5,k6,atmp,nlsolver,tab)
end

@cache mutable struct Kvaerno4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = Kvaerno4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  Kvaerno4ConstantCache(nlsolver,tab)
end

@cache mutable struct Kvaerno4Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
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

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = Kvaerno4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u); z₅ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Kvaerno4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,atmp,nlsolver,tab)
end

@cache mutable struct KenCarp4ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = KenCarp4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  KenCarp4ConstantCache(nlsolver,tab)
end

@cache mutable struct KenCarp4Cache{uType,rateType,uNoUnitsType,N,Tab,kType} <: SDIRKMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  z₅::uType
  z₆::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  k5::kType
  k6::kType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = KenCarp4Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  if typeof(f) <: SplitFunction
    k1 = similar(u); k2 = similar(u)
    k3 = similar(u); k4 = similar(u)
    k5 = similar(u); k6 = similar(u)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    k5 = nothing; k6 = nothing
    uf = UJacobianWrapper(f,t,p)
  end

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u); z₅ = zero(u)
  z₆ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  KenCarp4Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,k1,k2,k3,k4,k5,k6,atmp,nlsolver,tab)
end

@cache mutable struct Kvaerno5ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = Kvaerno5Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  Kvaerno5ConstantCache(nlsolver,tab)
end

@cache mutable struct Kvaerno5Cache{uType,rateType,uNoUnitsType,N,Tab} <: SDIRKMutableCache
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

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = Kvaerno5Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u); z₅ = zero(u)
  z₆ = zero(u); z₇ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  Kvaerno5Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,z₇,atmp,nlsolver,tab)
end

@cache mutable struct KenCarp5ConstantCache{N,Tab} <: OrdinaryDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  tab = KenCarp5Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))

  KenCarp5ConstantCache(nlsolver,tab)
end

@cache mutable struct KenCarp5Cache{uType,rateType,uNoUnitsType,N,Tab,kType} <: SDIRKMutableCache
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
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  k5::kType
  k6::kType
  k7::kType
  k8::kType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tab = KenCarp5Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ, tab.c3
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  if typeof(f) <: SplitFunction
    k1 = similar(u); k2 = similar(u)
    k3 = similar(u); k4 = similar(u)
    k5 = similar(u); k6 = similar(u)
    k7 = similar(u); k8 = similar(u)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    k5 = nothing; k6 = nothing
    k7 = nothing; k8 = nothing
  end

  z₁ = zero(u); z₂ = zero(u); z₃ = zero(u); z₄ = zero(u)
  z₅ = zero(u); z₆ = zero(u); z₇ = zero(u); z₈ = nlsolver.z
  atmp = similar(u,uEltypeNoUnits)

  KenCarp5Cache(u,uprev,fsalfirst,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,
                k1,k2,k3,k4,k5,k6,k7,k8,atmp,nlsolver,tab)
end
