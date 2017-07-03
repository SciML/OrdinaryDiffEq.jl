type ImplicitEulerCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  u_old::uArrayType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.u_old)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)
vecu_cache(c::ImplicitEulerCache) = (c.uhold,)
dual_cache(c::ImplicitEulerCache) = (c.dual_cache,)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  u_old = vec(tmp); k = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IE(f,u_old,t,t,dual_cache,size(u),eachindex(u))
  fsalfirst = zeros(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)

  ImplicitEulerCache{typeof(u),typeof(u_old),typeof(uhold),typeof(dual_cache),
                     typeof(k),typeof(rhs),typeof(nl_rhs)}(
                     u,uprev,uprev2,uhold,dual_cache,u_old,tmp,k,fsalfirst,rhs,nl_rhs)
end

immutable ImplicitEulerConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_IE_Scalar(f,u_old,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  ImplicitEulerConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,u_old,rhs,nl_rhs)
end

type TrapezoidCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  u_old::uArrayType
  fsalfirst::rateType
  dual_cache::DiffCacheType
  tmp::uType
  k::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.u_old)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)
vecu_cache(c::TrapezoidCache) = (c.uhold,)
dual_cache(c::TrapezoidCache) = (c.dual_cache,)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  u_old = vec(tmp); k = zeros(rate_prototype)
  uhold = vec(u); fsalfirst = zeros(rate_prototype)
  f_old = vec(fsalfirst)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  rhs = RHS_Trap(f,u_old,f_old,t,t,size(u),dual_cache,eachindex(u))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  TrapezoidCache{typeof(u),typeof(u_old),typeof(uhold),typeof(dual_cache),typeof(k),
    typeof(rhs),typeof(nl_rhs)}(u,uprev,uprev2,uhold,u_old,fsalfirst,dual_cache,tmp,k,rhs,nl_rhs)
end


immutable TrapezoidConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_Trap_Scalar(f,u_old,rate_prototype,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  TrapezoidConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,u_old,rhs,nl_rhs)
end

get_chunksize{uType,DiffCacheType,rateType,CS}(cache::ImplicitEulerCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::TrapezoidCache{uType,DiffCacheType,rateType,CS}) = CS
