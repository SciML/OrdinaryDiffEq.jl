struct GenericIIF1ConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

struct GenericIIF1Cache{uType,DiffCacheType,rhsType,nl_rhsType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  fsalfirst::rateType
  expA::expType
  k::rateType
end

u_cache(c::GenericIIF1Cache)    = ()
du_cache(c::GenericIIF1Cache)   = (c.rtmp1,c.fsalfirst,c.k)
dual_cache(c::GenericIIF1Cache) = (c.dual_cache,)

function alg_cache(alg::GenericIIF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  rhs = RHS_IIF_Scalar(f,zero(u),t,t,one(uEltypeNoUnits))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericIIF1ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::GenericIIF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  A = f.f1
  expA = expm(A*dt)
  rhs = RHS_IIF(f,tmp,t,t,uEltypeNoUnits(1//1),dual_cache)
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,u)
  GenericIIF1Cache(u,uprev,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,expA,k)
end

struct GenericIIF2ConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

struct GenericIIF2Cache{uType,DiffCacheType,rhsType,nl_rhsType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  fsalfirst::rateType
  expA::expType
  k::rateType
end

u_cache(c::GenericIIF2Cache)    = ()
du_cache(c::GenericIIF2Cache)   = (c.rtmp1,c.fsalfirst,c.k)
dual_cache(c::GenericIIF2Cache) = (c.dual_cache,)

function alg_cache(alg::GenericIIF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(u)
  rhs = RHS_IIF_Scalar(f,tmp,t,t,uEltypeNoUnits(1//2))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericIIF2ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::GenericIIF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  A = f.f1
  expA = expm(A*dt)
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  rhs = RHS_IIF(f,tmp,t,t,uEltypeNoUnits(1//2),dual_cache)
  nl_rhs = alg.nlsolve(Val{:init},rhs,u)
  GenericIIF2Cache(u,uprev,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,expA,k)
end

struct LawsonEulerCache{uType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  rtmp::rateType
  expA::expType
  fsalfirst::rateType
end

function alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  A = f.f1
  expA = expm(A*dt)
  LawsonEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),expA,zeros(rate_prototype))
end

u_cache(c::LawsonEulerCache) = ()
du_cache(c::LawsonEulerCache) = (c.k,c.fsalfirst,c.rtmp)

struct LawsonEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = LawsonEulerConstantCache()

struct NorsettEulerCache{uType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  rtmp::rateType
  expA::expType
  phi1::expType
  fsalfirst::rateType
end

function alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  A = f.f1
  expA = expm(A*dt)
  phi1 = ((expA-I)/A)
  NorsettEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),expA,phi1,zeros(rate_prototype))
end

u_cache(c::NorsettEulerCache) = ()
du_cache(c::NorsettEulerCache) = (c.k,c.fsalfirst)

struct NorsettEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = NorsettEulerConstantCache()
