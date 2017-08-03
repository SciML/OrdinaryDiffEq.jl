type ImplicitEulerCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  C::uArrayType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.C)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)
vecu_cache(c::ImplicitEulerCache) = (c.uhold,)
dual_cache(c::ImplicitEulerCache) = (c.dual_cache,)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  C = vec(tmp); k = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = ImplicitRHS(f,C,t,t,t,dual_cache)
  fsalfirst = zeros(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)

  ImplicitEulerCache{typeof(u),typeof(C),typeof(uhold),typeof(dual_cache),
                     typeof(k),typeof(rhs),typeof(nl_rhs)}(
                     u,uprev,uprev2,uhold,dual_cache,C,tmp,k,fsalfirst,rhs,nl_rhs)
end

struct ImplicitEulerConstantCache{F} <: OrdinaryDiffEqConstantCache
  uf::F
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ImplicitEulerConstantCache(uf)
end

type TrapezoidCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  C::uArrayType
  fsalfirst::rateType
  dual_cache::DiffCacheType
  tmp::uType
  k::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.C)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)
vecu_cache(c::TrapezoidCache) = (c.uhold,)
dual_cache(c::TrapezoidCache) = (c.dual_cache,)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  C = vec(tmp); k = zeros(rate_prototype)
  uhold = vec(u); fsalfirst = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  rhs = ImplicitRHS(f,C,t,t,t,dual_cache)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  TrapezoidCache{typeof(u),typeof(C),typeof(uhold),typeof(dual_cache),typeof(k),
    typeof(rhs),typeof(nl_rhs)}(u,uprev,uprev2,uhold,C,fsalfirst,dual_cache,tmp,k,rhs,nl_rhs)
end


struct TrapezoidConstantCache{F} <: OrdinaryDiffEqConstantCache
  uf::F
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  TrapezoidConstantCache(uf)
end
