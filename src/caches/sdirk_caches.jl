mutable struct ImplicitEulerCache{uType,rateType,J,JC,UF,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  tmp::uType
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
  dz = similar(u); tmp = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,tmp,J,W,jac_config,uf,ηold,κ,tol,10000)
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
  tmp::uType
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
  dz = similar(u); tmp = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,tmp,J,W,jac_config,uf,ηold,κ,tol,10000,uprev3,tprev2)
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
  tmp::uType
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
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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
              u,uprev,uᵧ,du1,fsalfirst,k,zprev,zᵧ,z,Δzᵧ,Δz,tmp,J,
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
  tmp::uType
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
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz₁,dz₂,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000)
end

mutable struct SSPSDIRK2ConstantCache{F,uEltypeNoUnits,uType,tType} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
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

  SSPSDIRK2ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2)
end

mutable struct SSPSDIRK2Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  dz₁::uType
  dz₂::uType
  tmp::uType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::SSPSDIRK2Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::SSPSDIRK2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SSPSDIRK2,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  SSPSDIRK2Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,dz₁,dz₂,tmp,J,
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
  tmp::uType
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
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz₁,dz₂,dz₃,dz₄,tmp,J,
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
  tmp::uType
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
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz₁,dz₂,dz₃,dz₄,tmp,J,
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
  z₅::uType
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  dz₅::uType
  tmp::uType
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
  z₅ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  dz₅ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz₁,dz₂,dz₃,dz₄,dz₅,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct Hairer4ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
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

  if typeof(alg) <: Hairer4
    tab = Hairer4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))
  else
    tab = Hairer42Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))
  end

  Hairer4ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct Hairer4Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
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
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  dz₅::uType
  tmp::uType
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

u_cache(c::Hairer4Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::Hairer4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Union{Hairer4,Hairer42},u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  dz₅ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  if typeof(alg) <: Hairer4
    tab = Hairer4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))
  else # Hairer42
    tab = Hairer42Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))
  end

  ηold = one(uEltypeNoUnits)

  Hairer4Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz₁,dz₂,dz₃,dz₄,dz₅,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno4ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
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

  tab = Kvaerno4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  Kvaerno4ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct Kvaerno4Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
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
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  dz₅::uType
  tmp::uType
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

u_cache(c::Kvaerno4Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::Kvaerno4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno4,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  dz₅ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  tab = Kvaerno4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  Kvaerno4Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,dz₁,dz₂,dz₃,dz₄,dz₅,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp4ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
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

  tab = KenCarp4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  KenCarp4ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct KenCarp4Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
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
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  dz₅::uType
  dz₆::uType
  tmp::uType
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

u_cache(c::KenCarp4Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::KenCarp4Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::KenCarp4,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u); z₆ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  dz₅ = similar(u); dz₆ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  tab = KenCarp4Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  KenCarp4Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,dz₁,dz₂,dz₃,dz₄,dz₅,dz₆,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct Kvaerno5ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
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

  tab = Kvaerno5Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  Kvaerno5ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct Kvaerno5Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
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
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  dz₅::uType
  dz₆::uType
  dz₇::uType
  tmp::uType
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

u_cache(c::Kvaerno5Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::Kvaerno5Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Kvaerno5,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u); z₆ = similar(u); z₇ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  dz₅ = similar(u); dz₆ = similar(u); dz₇ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  tab = Kvaerno5Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  Kvaerno5Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,dz₁,dz₂,dz₃,dz₄,dz₅,dz₆,dz₇,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end

mutable struct KenCarp5ConstantCache{F,uEltypeNoUnits,uType,tType,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  uprev3::uType
  tprev2::tType
  tab::Tab
end

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
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

  tab = KenCarp5Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  KenCarp5ConstantCache(uf,ηold,κ,tol,10000,uprev3,tprev2,tab)
end

mutable struct KenCarp5Cache{uType,rateType,J,JC,UF,uEltypeNoUnits,tType,Tab} <: OrdinaryDiffEqMutableCache
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
  dz₁::uType
  dz₂::uType
  dz₃::uType
  dz₄::uType
  dz₅::uType
  dz₆::uType
  dz₇::uType
  dz₈::uType
  tmp::uType
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

u_cache(c::KenCarp5Cache)    = (c.uprev2,c.zᵧ,c.z,c.Δzᵧ,c.Δz)
du_cache(c::KenCarp5Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::KenCarp5,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u); z₂ = similar(u); z₃ = similar(u); z₄ = similar(u)
  z₅ = similar(u); z₆ = similar(u); z₇ = similar(u); z₈ = similar(u)
  dz₁ = similar(u); dz₂ = similar(u); dz₃ = similar(u); dz₄ = similar(u)
  dz₅ = similar(u); dz₆ = similar(u); dz₇ = similar(u); dz₈ = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype); tmp = similar(u)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t,vec(uprev),vec(du1))
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

  tab = KenCarp5Tableau(real(uEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  KenCarp5Cache{typeof(u),typeof(rate_prototype),typeof(J),typeof(jac_config),
              typeof(uf),uEltypeNoUnits,typeof(t),typeof(tab)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,z₅,z₆,z₇,z₈,
              dz₁,dz₂,dz₃,dz₄,dz₅,dz₆,dz₇,dz₈,tmp,J,
              W,jac_config,uf,ηold,κ,tol,10000,tab)
end
