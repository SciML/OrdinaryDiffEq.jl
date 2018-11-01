@cache mutable struct KenCarp3ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  tab = KenCarp3Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  KenCarp3ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct KenCarp3Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F,kType} <: SDIRKMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
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

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u)
  z₃ = similar(u); z₄ = z
  atmp = similar(u,uEltypeNoUnits)

  if typeof(f) <: SplitFunction
    k1 = similar(u); k2 = similar(u)
    k3 = similar(u); k4 = similar(u)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end
  tab = KenCarp3Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  KenCarp3Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,k1,k2,k3,k4,dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolve,tab)
end

@cache mutable struct Kvaerno4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  uprev3 = u
  tprev2 = t

  tab = Kvaerno4Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  Kvaerno4ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct Kvaerno4Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F} <: SDIRKMutableCache
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

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = similar(u)
  z₅ = z
  atmp = similar(u,uEltypeNoUnits)
  tab = Kvaerno4Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  Kvaerno4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolve,tab)
end

@cache mutable struct KenCarp4ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  uprev3 = u
  tprev2 = t
  tab = KenCarp4Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  KenCarp4ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct KenCarp4Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F,kType} <: SDIRKMutableCache
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
  z₆::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  k5::kType
  k6::kType
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

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u)
  z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u); z₆ = z
  atmp = similar(u,uEltypeNoUnits)

  if typeof(f) <: SplitFunction
    k1 = similar(u); k2 = similar(u)
    k3 = similar(u); k4 = similar(u)
    k5 = similar(u); k6 = similar(u)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    k5 = nothing; k6 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end

  tab = KenCarp4Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  KenCarp4Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,k1,k2,k3,k4,k5,k6,
                dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolve,tab)
end

@cache mutable struct Kvaerno5ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  tab = Kvaerno5Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  Kvaerno5ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct Kvaerno5Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F} <: SDIRKMutableCache
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
  z₆::uType
  z₇::uType
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

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u); z₆ = similar(u);
  z₇ = z
  atmp = similar(u,uEltypeNoUnits)
  tab = Kvaerno5Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  Kvaerno5Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolve,tab)
end

@cache mutable struct KenCarp5ConstantCache{F,N,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  @oopnlcachefields
  tab = KenCarp5Tableau(uToltype,real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  KenCarp5ConstantCache(uf,nlsolve,tab)
end

@cache mutable struct KenCarp5Cache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F,kType} <: SDIRKMutableCache
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

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  @iipnlcachefields
  z₁ = similar(u); z₂ = similar(u);
  z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u); z₆ = similar(u);
  z₇ = similar(u); z₈ = z
  atmp = similar(u,uEltypeNoUnits)

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
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end

  tab = KenCarp5Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))

  KenCarp5Cache(u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,
                k1,k2,k3,k4,k5,k6,k7,k8,
                dz,b,tmp,atmp,J,
                W,uf,jac_config,linsolve,nlsolve,tab)
end
