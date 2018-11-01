@cache mutable struct GenericImplicitEulerCache{uType,DiffCacheType,uNoUnitsType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  dual_cache::DiffCacheType
  tmp::uType
  atmp::uNoUnitsType
  k::rateType
  fsalfirst::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::GenericImplicitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  k = zero(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  rhs = ImplicitRHS(f,tmp,t,t,t,dual_cache,p)
  fsalfirst = zero(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,u)

  GenericImplicitEulerCache{typeof(u),typeof(dual_cache),
                     typeof(atmp),typeof(k),typeof(rhs),typeof(nl_rhs)}(
                     u,uprev,uprev2,dual_cache,tmp,atmp,k,fsalfirst,rhs,nl_rhs)
end

struct GenericImplicitEulerConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::GenericImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uBottomEltypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uhold = Vector{typeof(u)}(undef, 1)
  rhs = ImplicitRHS_Scalar(f,zero(u),t,t,t,p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericImplicitEulerConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,rhs,nl_rhs)
end

@cache mutable struct GenericTrapezoidCache{uType,DiffCacheType,uNoUnitsType,rateType,rhsType,nl_rhsType,tType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  dual_cache::DiffCacheType
  tmp::uType
  atmp::uNoUnitsType
  k::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::GenericTrapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  rhs = ImplicitRHS(f,tmp,t,t,t,dual_cache,p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,u)
  uprev3 = similar(u)
  tprev2 = t
  GenericTrapezoidCache{typeof(u),typeof(dual_cache),typeof(atmp),typeof(k),
                        typeof(rhs),typeof(nl_rhs),typeof(t)}(
                        u,uprev,uprev2,fsalfirst,
                        dual_cache,tmp,atmp,k,rhs,nl_rhs,uprev3,tprev2)
end


@cache mutable struct GenericTrapezoidConstantCache{vecuType,rhsType,nl_rhsType,uType,tType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::GenericTrapezoid,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uhold = Vector{typeof(u)}(undef, 1)
  rhs = ImplicitRHS_Scalar(f,zero(u),t,t,t,p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  uprev3 = u
  tprev2 = t
  GenericTrapezoidConstantCache(uhold,rhs,nl_rhs,uprev3,tprev2)
end

get_chunksize(cache::GenericImplicitEulerCache{uType,DiffCacheType,rateType,CS}) where {uType,DiffCacheType,rateType,CS} = CS
get_chunksize(cache::GenericTrapezoidCache{uType,DiffCacheType,rateType,CS}) where {uType,DiffCacheType,rateType,CS} = CS
