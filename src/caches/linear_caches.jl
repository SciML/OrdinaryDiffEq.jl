struct LinearImplicitEulerCache{uType,rateType,J} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  W::J
  k::rateType
end

u_cache(c::LinearImplicitEulerCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::LinearImplicitEulerCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::LinearImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  W = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  k = zeros(rate_prototype); fsalfirst = zeros(rate_prototype)
  LinearImplicitEulerCache(u,uprev,uprev2,fsalfirst,W,k)
end

struct LinearImplicitEulerConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(alg::LinearImplicitEuler,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  LinearImplicitEulerConstantCache()
end

struct StrangSplittingCache{uType,rateType,J} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  W::J
  k::rateType
end

u_cache(c::StrangSplittingCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::StrangSplittingCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::StrangSplitting,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  W = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  k = zeros(rate_prototype); fsalfirst = zeros(rate_prototype)
  StrangSplittingCache(u,uprev,uprev2,fsalfirst,W,k)
end

struct StrangSplittingConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(alg::StrangSplitting,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  StrangSplittingConstantCache()
end
