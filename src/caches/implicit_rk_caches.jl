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

u_cache(c::TRBDF2Cache)    = (c.uprev2,c.z,c.dz)
du_cache(c::TRBDF2Cache)   = (c.k,c.fsalfirst)

function alg_cache(alg::TRBDF2,u,rate_prototype,uEltypeNoUnits,
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

  TRBDF2Cache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,J,W,jac_config,uf,ηold,κ,tol,10000,uprev3,tprev2)
end
