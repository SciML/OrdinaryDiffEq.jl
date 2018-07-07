mutable struct KenCarp3ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  if typeof(f) <: SplitFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
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

  tab = KenCarp3Tableau(uToltype,real(tTypeNoUnits))

  KenCarp3ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp3Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F,kType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::KenCarp3Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.dz)
du_cache(c::KenCarp3Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u))
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
  atmp = similar(u,uEltypeNoUnits,axes(u))

  if typeof(f) <: SplitFunction
    k1 = similar(u,axes(u)); k2 = similar(u,axes(u))
    k3 = similar(u,axes(u)); k4 = similar(u,axes(u))
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end
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

  tab = KenCarp3Tableau(uToltype,real(tTypeNoUnits))

  ηold = one(uToltype)

  KenCarp3Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve),typeof(k1)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,k1,k2,k3,k4,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno4ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  tab = Kvaerno4Tableau(uToltype,real(tTypeNoUnits))

  Kvaerno4ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno4Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F} <: OrdinaryDiffEqMutableCache
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

u_cache(c::Kvaerno4Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.dz₁,c.dz₂,c.dz₃,c.dz₄,c.dz₅)
du_cache(c::Kvaerno4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
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

  tab = Kvaerno4Tableau(uToltype,real(tTypeNoUnits))

  ηold = one(uToltype)

  Kvaerno4Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp4ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  if typeof(f) <: SplitFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
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

  tab = KenCarp4Tableau(uToltype,real(tTypeNoUnits))

  KenCarp4ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp4Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F,kType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::KenCarp4Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.z₆,c.dz)
du_cache(c::KenCarp4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u))
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  z₅ = similar(u,axes(u)); z₆ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
  atmp = similar(u,uEltypeNoUnits,axes(u))

  if typeof(f) <: SplitFunction
    k1 = similar(u,axes(u)); k2 = similar(u,axes(u))
    k3 = similar(u,axes(u)); k4 = similar(u,axes(u))
    k5 = similar(u,axes(u)); k6 = similar(u,axes(u))
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    k5 = nothing; k6 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end
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

  tab = KenCarp4Tableau(uToltype,real(tTypeNoUnits))

  ηold = one(uToltype)

  KenCarp4Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve),typeof(k1)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,k1,k2,k3,k4,k5,k6,
              dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno5ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
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

  tab = Kvaerno5Tableau(uToltype,real(tTypeNoUnits))

  Kvaerno5ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno5Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F} <: OrdinaryDiffEqMutableCache
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

u_cache(c::Kvaerno5Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.z₆,c.z₇,c.dz)
du_cache(c::Kvaerno5Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u));
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  z₅ = similar(u,axes(u)); z₆ = similar(u,axes(u));
  z₇ = similar(u,axes(u))
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

  tab = Kvaerno5Tableau(uToltype,real(tTypeNoUnits))

  ηold = one(uToltype)

  Kvaerno5Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp5ConstantCache{F,uToltype,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uToltype
  κ::uToltype
  tol::uToltype
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  if typeof(f) <: SplitFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
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

  tab = KenCarp5Tableau(uToltype,real(tTypeNoUnits))

  KenCarp5ConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp5Cache{uType,rateType,uNoUnitsType,J,UF,JC,uToltype,Tab,F,kType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::KenCarp5Cache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.z₅,c.z₆,c.z₇,c.z₈,c.dz)
du_cache(c::KenCarp5Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zero(rate_prototype)
  J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u));
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  z₅ = similar(u,axes(u)); z₆ = similar(u,axes(u));
  z₇ = similar(u,axes(u)); z₈ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = similar(u); b = similar(u,axes(u));
  atmp = similar(u,uEltypeNoUnits,axes(u))

  if typeof(f) <: SplitFunction
    k1 = similar(u,axes(u)); k2 = similar(u,axes(u))
    k3 = similar(u,axes(u)); k4 = similar(u,axes(u))
    k5 = similar(u,axes(u)); k6 = similar(u,axes(u))
    k7 = similar(u,axes(u)); k8 = similar(u,axes(u))
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    k5 = nothing; k6 = nothing
    k7 = nothing; k8 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end

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

  tab = KenCarp5Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uToltype)

  KenCarp5Cache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uToltype,typeof(tab),typeof(linsolve),typeof(k1)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,
              k1,k2,k3,k4,k5,k6,k7,k8,
              dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end
