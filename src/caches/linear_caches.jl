struct LinearImplicitEulerCache{uType,rateType,J,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  W::J
  k::rateType
  tmp::uType
  atmp::uNoUnitsType
end

u_cache(c::LinearImplicitEulerCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::LinearImplicitEulerCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::LinearImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,::Type{Val{true}})
  W = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  k = zeros(rate_prototype); fsalfirst = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u))
  LinearImplicitEulerCache(u,uprev,uprev2,fsalfirst,W,k,tmp,atmp)
end

struct LinearImplicitEulerConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(alg::LinearImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,::Type{Val{false}})
  LinearImplicitEulerConstantCache()
end

struct MidpointSplittingCache{uType,rateType,J} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  tmp::uType
  fsalfirst::rateType
  W::J
  k::rateType
end

u_cache(c::MidpointSplittingCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::MidpointSplittingCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::MidpointSplitting,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,::Type{Val{true}})
  W = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  k = zeros(rate_prototype); fsalfirst = zeros(rate_prototype)
  MidpointSplittingCache(u,uprev,uprev2,similar(u),fsalfirst,W,k)
end

struct MidpointSplittingConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(alg::MidpointSplitting,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,::Type{Val{false}})
  MidpointSplittingConstantCache()
end
