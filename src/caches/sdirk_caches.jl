mutable struct ImplicitEulerCache{uType,rateType,J,JC,UF,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u)
  dz = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,J,W,jac_config,uf,ηold,κ,tol,10000)
end

mutable struct ImplicitEulerConstantCache{F,uEltypeNoUnits} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  ImplicitEulerConstantCache(uf,ηold,κ,tol,100000)
end

mutable struct TrapezoidConstantCache{F,uEltypeNoUnits,uType,tType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  TrapezoidConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct TrapezoidCache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u)
  dz = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  uprev3 = similar(u)
  tprev2 = t

  ηold = one(uEltypeNoUnits)

  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,J,W,jac_config,uf,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct TRBDF2ConstantCache{F,uEltypeNoUnits,uType,tType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  TRBDF2ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct TRBDF2Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uᵧ::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  zprev::uType
  zᵧ::uType
  z::uType
  Δzᵧ::uType
  Δz::uType
  est::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::TRBDF2Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::TRBDF2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  uᵧ = similar(u)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  zprev = similar(u); zᵧ = similar(u); z = similar(u)
  Δzᵧ = similar(u); Δz = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); est = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  ηold = one(uEltypeNoUnits)

  TRBDF2Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t)}(
              u,uprev,uᵧ,du1,fsalfirst,k,zprev,zᵧ,z,Δzᵧ,Δz,est,J,
              W,jac_config,uf,ηold,κ,tol,10000)
end

mutable struct SDIRK2ConstantCache{F,uEltypeNoUnits,uType,tType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  SDIRK2ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct SDIRK2Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  dz₁::uType
  dz₂::uType
  est::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::SDIRK2Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::SDIRK2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SDIRK2,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); est = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  ηold = one(uEltypeNoUnits)

  SDIRK2Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz₁,dz₂,est,J,
              W,jac_config,uf,ηold,κ,tol,10000)
end

mutable struct Kvaerno3ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = Kvaerno3Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  Kvaerno3ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct Kvaerno3Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  est::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
end

u_cache(c::Kvaerno3Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::Kvaerno3Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno3,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); est = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = Kvaerno3Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  Kvaerno3Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz₁,dz₂,dz₃,dz₄,est,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp3ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = KenCarp3Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  KenCarp3ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct KenCarp3Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  est::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
end

u_cache(c::KenCarp3Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::KenCarp3Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::KenCarp3,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); est = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = KenCarp3Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  KenCarp3Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz₁,dz₂,dz₃,dz₄,est,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct Cash4ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = Cash4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  Cash4ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct Cash4Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  est::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
end

u_cache(c::Cash4Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::Cash4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Cash4,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); est = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
  end

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = Cash4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  Cash4Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz₁,dz₂,dz₃,dz₄,est,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end
