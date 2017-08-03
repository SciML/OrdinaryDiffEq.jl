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
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})

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
  ImplicitEulerCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,J,W,jac_config,uf,ηold)
end

mutable struct ImplicitEulerConstantCache{F,uEltypeNoUnits} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  ImplicitEulerConstantCache(uf,ηold)
end

mutable struct TrapezoidConstantCache{F,uEltypeNoUnits} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)
  TrapezoidConstantCache(uf,ηold)
end

mutable struct TrapezoidCache{uType,rateType,J,JC,UF,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
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
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})

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

  TrapezoidCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,J,W,jac_config,uf,ηold)
end
