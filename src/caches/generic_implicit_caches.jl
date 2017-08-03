mutable struct GenericImplicitEulerCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::GenericImplicitEulerCache)    = (c.uprev2,c.C)
du_cache(c::GenericImplicitEulerCache)   = (c.k,c.fsalfirst)
vecu_cache(c::GenericImplicitEulerCache) = (c.uhold,)
dual_cache(c::GenericImplicitEulerCache) = (c.dual_cache,)

function alg_cache(alg::GenericImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  tmp = similar(u)
  C = vec(tmp); k = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = ImplicitRHS(f,C,t,t,t,dual_cache)
  fsalfirst = zeros(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)

  GenericImplicitEulerCache{typeof(u),typeof(C),typeof(uhold),typeof(dual_cache),
                     typeof(k),typeof(rhs),typeof(nl_rhs)}(
                     u,uprev,uprev2,uhold,dual_cache,C,tmp,k,fsalfirst,rhs,nl_rhs)
end

struct GenericImplicitEulerConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  C::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::GenericImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  C = Vector{typeof(u)}(1)
  rhs = ImplicitRHS_Scalar(f,C,t,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericImplicitEulerConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,C,rhs,nl_rhs)
end

mutable struct GenericTrapezoidCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
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

u_cache(c::GenericTrapezoidCache)    = (c.uprev2,c.C)
du_cache(c::GenericTrapezoidCache)   = (c.k,c.fsalfirst)
vecu_cache(c::GenericTrapezoidCache) = (c.uhold,)
dual_cache(c::GenericTrapezoidCache) = (c.dual_cache,)

function alg_cache(alg::GenericTrapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  tmp = similar(u)
  C = vec(tmp); k = zeros(rate_prototype)
  uhold = vec(u); fsalfirst = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  rhs = ImplicitRHS(f,C,t,t,t,dual_cache)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericTrapezoidCache{typeof(u),typeof(C),typeof(uhold),typeof(dual_cache),typeof(k),
    typeof(rhs),typeof(nl_rhs)}(u,uprev,uprev2,uhold,C,fsalfirst,dual_cache,tmp,k,rhs,nl_rhs)
end


struct GenericTrapezoidConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  C::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::GenericTrapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  C = Vector{typeof(u)}(1)
  rhs = ImplicitRHS_Scalar(f,C,t,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericTrapezoidConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,C,rhs,nl_rhs)
end

get_chunksize{uType,DiffCacheType,rateType,CS}(cache::GenericImplicitEulerCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::GenericTrapezoidCache{uType,DiffCacheType,rateType,CS}) = CS
