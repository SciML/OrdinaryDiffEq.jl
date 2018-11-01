@cache struct MidpointSplittingCache{uType,rateType,WType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  tmp::uType
  fsalfirst::rateType
  W::WType
  k::rateType
end

function alg_cache(alg::MidpointSplitting,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  W = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
  k = zero(rate_prototype); fsalfirst = zero(rate_prototype)
  MidpointSplittingCache(u,uprev,uprev2,similar(u),fsalfirst,W,k)
end

struct MidpointSplittingConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(alg::MidpointSplitting,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  MidpointSplittingConstantCache()
end

struct LinearExponentialConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::LinearExponential,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = LinearExponentialConstantCache()

@cache struct LinearExponentialCache{uType,rateType,KsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  KsCache::KsType # different depending on alg.krylov
end

function alg_cache(alg::LinearExponential,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  rtmp = zero(rate_prototype)
  n = length(u)
  T = eltype(u)
  m = min(alg.m, n)

  if alg.krylov == :off
    KsCache = nothing
  elseif alg.krylov == :simple
    Ks = KrylovSubspace{T}(n, m)
    expv_cache = ExpvCache{T}(m)
    KsCache = (Ks, expv_cache)
  elseif alg.krylov == :adaptive
    KsCache = _phiv_timestep_caches(u, m, 0)
  else
    throw(ArgumentError("Unknown krylov setting $(alg.krylov). Can be :off, :simple or :adaptive."))
  end
  LinearExponentialCache(u,uprev,tmp,rtmp,KsCache)
end
