struct IIF1ConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

struct IIF1Cache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  fsalfirst::rateType
  expA::expType
  k::rateType
end

u_cache(c::IIF1Cache)    = (c.uprev2,c.u_old)
du_cache(c::IIF1Cache)   = (c.rtmp1,c.tmp,c.fsalfirst,c.k)
vecu_cache(c::IIF1Cache) = (c.uhold,)
dual_cache(c::IIF1Cache) = (c.dual_cache,)

function alg_cache(alg::IIF1,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(rate_prototype)
  rhs = RHS_IIF_Scalar(f,tmp,t,t,one(uEltypeNoUnits))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF1ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF1,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  A = f.f1
  expA = expm(A*dt)
  rhs = RHS_IIF(f,tmp,t,t,dual_cache,uEltypeNoUnits(1//1))
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF1Cache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,expA,k)
end

struct IIF2ConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

struct IIF2Cache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  fsalfirst::rateType
  expA::expType
  k::rateType
end

u_cache(c::IIF2Cache)    = (c.uprev2,c.u_old)
du_cache(c::IIF2Cache)   = (c.rtmp1,c.tmp,c.fsalfirst,c.k)
vecu_cache(c::IIF2Cache) = (c.uhold,)
dual_cache(c::IIF2Cache) = (c.dual_cache,)

function alg_cache(alg::IIF2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(rate_prototype)
  rhs = RHS_IIF_Scalar(f,tmp,t,t,uEltypeNoUnits(1//2))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF2ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  A = f.f1
  expA = expm(A*dt)
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  rhs = RHS_IIF(f,tmp,t,t,dual_cache,uEltypeNoUnits(1//2))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF2Cache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,expA,k)
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

function alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  A = f.f1
  expA = expm(A*dt)
  LawsonEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),expA,zeros(rate_prototype))
end

u_cache(c::LawsonEulerCache) = ()
du_cache(c::LawsonEulerCache) = (c.k,c.fsalfirst,c.rtmp)

struct LawsonEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = LawsonEulerConstantCache()

struct NorsettEulerCache{uType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  rtmp::rateType
  expA::expType
  fsalfirst::rateType
end

function alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  A = f.f1
  expA = expm(A*dt)
  NorsettEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),expA,zeros(rate_prototype))
end

u_cache(c::NorsettEulerCache) = ()
du_cache(c::NorsettEulerCache) = (c.k,c.fsalfirst)

struct NorsettEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = NorsettEulerConstantCache()
