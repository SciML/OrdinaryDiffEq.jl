######################################
# IIF Caches

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

function alg_cache(alg::GenericIIF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  rhs = RHS_IIF_Scalar(f,zero(u),t,t,one(uEltypeNoUnits),p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericIIF1ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::GenericIIF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  A = f.f1
  expA = expm(A*dt)
  rhs = RHS_IIF(f,tmp,t,t,uEltypeNoUnits(1//1),dual_cache,p)
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

function alg_cache(alg::GenericIIF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(u)
  rhs = RHS_IIF_Scalar(f,tmp,t,t,uEltypeNoUnits(1//2),p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericIIF2ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::GenericIIF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  A = f.f1
  expA = expm(A*dt)
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  rhs = RHS_IIF(f,tmp,t,t,uEltypeNoUnits(1//2),dual_cache,p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,u)
  GenericIIF2Cache(u,uprev,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,expA,k)
end

#################################################
# Classical ExpRK method caches
abstract type ExpRKCache <: OrdinaryDiffEqMutableCache end
abstract type ExpRKConstantCache <: OrdinaryDiffEqConstantCache end

# Precomputation of exponential-like operators
expRK_operators(::LawsonEuler, dt, A) = expm(dt * A)
expRK_operators(::NorsettEuler, dt, A) = phi(dt * A, 1)[2]
function expRK_operators(::ETDRK2, dt, A)
  P = phi(dt * A, 2)
  return P[2], P[3] # ϕ1(hA), ϕ2(hA)
end
function expRK_operators(::ETDRK3, dt, A)
  Phalf = phi(dt/2 * A, 1)
  A21 = 0.5Phalf[2]
  P = phi(dt * A, 3)
  phi1, phi2, phi3 = P[2], P[3], P[4]
  A3 = phi1
  B1 = 4phi3 - 3phi2 + phi1
  B2 = -8phi3 + 4phi2
  B3 = 4phi3 - phi2
  return A21, A3, B1, B2, B3
end
function expRK_operators(::ETDRK4, dt, A)
  P = phi(dt * A, 3)
  Phalf = phi(dt/2 * A, 1)
  A21 = 0.5Phalf[2] # A32 = A21
  A41 = (dt/4 * A) * Phalf[2]^2
  A43 = Phalf[2]
  B1 = P[2] - 3P[3] + 4P[4]
  B2 = 2P[3] - 4P[4] # B3 = B2
  B4 = -P[3] + 4P[4]
  return A21, A41, A43, B1, B2, B4
end

# Unified constructor for constant caches
for (Alg, Cache) in [(:LawsonEuler, :LawsonEulerConstantCache), 
                     (:NorsettEuler, :NorsettEulerConstantCache),
                     (:ETDRK2, :ETDRK2ConstantCache),
                     (:ETDRK3, :ETDRK3ConstantCache),
                     (:ETDRK4, :ETDRK4ConstantCache)]
  @eval struct $Cache{opType} <: ExpRKConstantCache
    ops::opType # precomputed operators
  end

  @eval function alg_cache(alg::$Alg,u,rate_prototype,uEltypeNoUnits,
    uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
    if alg.krylov
      ops = nothing # no caching
    else
      isa(f, SplitFunction) || throw(ArgumentError("Caching can only be used with SplitFunction"))
      A = isa(f.f1, DiffEqArrayOperator) ? f.f1.A * f.f1.α.coeff : full(f.f1)
      ops = expRK_operators(alg, dt, A)
    end
    return $Cache(ops)
  end
end

struct LawsonEulerCache{uType,rateType,JType,expType,KsType,KsCacheType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  G::rateType
  Jcache::JType
  exphA::expType
  Ks::KsType
  KsCache::KsCacheType
end

function alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if isa(f, SplitFunction)
    Jcache = nothing
  else
    # TODO: sparse Jacobian support
    Jcache = Matrix{eltype(u)}(length(u), length(u))
  end
  if alg.krylov
    exphA = nothing # no caching
    m = min(alg.m, length(u))
    Ks = KrylovSubspace{eltype(u)}(length(u), m)
    KsCache = Matrix{eltype(u)}(m, m)
  else
    Ks = nothing
    KsCache = nothing
    A = f.f1
    if isa(A, DiffEqArrayOperator)
      _A = A.A * A.α.coeff
    else
      _A = full(A)
    end
    exphA = expRK_operators(alg, dt, _A)
  end
  LawsonEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),Jcache,exphA,Ks,KsCache)
end

u_cache(c::LawsonEulerCache) = ()
du_cache(c::LawsonEulerCache) = (c.rtmp)

struct NorsettEulerCache{uType,rateType,JType,expType,KsType,KsCacheType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  G::rateType
  Jcache::JType
  phihA::expType
  Ks::KsType
  KsCache::KsCacheType
end

function alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if isa(f, SplitFunction)
    Jcache = nothing
  else
    # TODO: sparse Jacobian support
    Jcache = Matrix{eltype(u)}(length(u), length(u))
  end
  if alg.krylov
    phihA = nothing # no caching
    n = length(u)
    m = min(alg.m, length(u))
    T = eltype(u)
    Ks = KrylovSubspace{T}(n, m)
    KsCache = (Matrix{T}(n, 2), Vector{T}(m), Matrix{T}(m, m), Matrix{T}(m + 1, m + 1), 
      Matrix{T}(m, 2))
  else
    Ks = nothing
    KsCache = nothing
    A = f.f1
    if isa(A, DiffEqArrayOperator)
      _A = A.A * A.α.coeff
    else
      _A = full(A)
    end
    phihA = expRK_operators(alg, dt, _A)
  end
  NorsettEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),Jcache,phihA,Ks,KsCache)
end

u_cache(c::NorsettEulerCache) = ()
du_cache(c::NorsettEulerCache) = (c.rtmp)

struct ETDRK2Cache{uType,rateType,JType,opType,KsType,KsCacheType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  F2::rateType
  Jcache::JType
  ops::opType
  Ks::KsType
  KsCache::KsCacheType
end

function alg_cache(alg::ETDRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if isa(f, SplitFunction)
    Jcache = nothing
  else
    # TODO: sparse Jacobian support
    Jcache = Matrix{eltype(u)}(length(u), length(u))
  end
  if alg.krylov
    ops = nothing # no caching
    n = length(u)
    m = min(alg.m, length(u))
    T = eltype(u)
    Ks = KrylovSubspace{T}(n, m)
    w1 = Matrix{T}(n, 3)
    w2 = Matrix{T}(n, 3)
    phiv_caches = (Vector{T}(m), Matrix{T}(m, m), Matrix{T}(m + 2, m + 2), Matrix{T}(m, 3))
    KsCache = (w1, w2, phiv_caches)
  else
    Ks = nothing
    KsCache = nothing
    A = f.f1
    if isa(A, DiffEqArrayOperator)
      _A = A.A * A.α.coeff
    else
      _A = full(A)
    end
    ops = expRK_operators(alg, dt, _A)
  end
  ETDRK2Cache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),Jcache,ops,Ks,KsCache)
end

struct ETDRK3Cache{uType,rateType,JType,opType,KsType,KsCacheType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  Au::rateType
  F2::rateType
  F3::rateType
  Jcache::JType
  ops::opType
  Ks::KsType
  KsCache::KsCacheType
end

function alg_cache(alg::ETDRK3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if isa(f, SplitFunction)
    Jcache = nothing
  else
    # TODO: sparse Jacobian support
    Jcache = Matrix{eltype(u)}(length(u), length(u))
  end
  if alg.krylov
    ops = nothing # no caching
    n = length(u)
    m = min(alg.m, length(u))
    T = eltype(u)
    Ks = KrylovSubspace{T}(n, m)
    w1_half = Matrix{T}(n, 2); w1 = Matrix{T}(n, 4); w2 = Matrix{T}(n, 4); w3 = Matrix{T}(n, 4)
    phiv_caches = (Vector{T}(m), Matrix{T}(m, m), Matrix{T}(m + 3, m + 3), Matrix{T}(m, 4))
    KsCache = (w1_half, w1, w2, w3, phiv_caches)
  else
    Ks = nothing
    KsCache = nothing
    A = f.f1
    if isa(A, DiffEqArrayOperator)
      _A = A.A * A.α.coeff
    else
      _A = full(A)
    end
    ops = expRK_operators(alg, dt, _A)
  end
  ETDRK3Cache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),
    zeros(rate_prototype),zeros(rate_prototype),Jcache,ops,Ks,KsCache)
end

struct ETDRK4Cache{uType,rateType,JType,opType,KsType,KsCacheType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  Au::rateType
  F2::rateType
  F3::rateType
  F4::rateType
  Jcache::JType
  ops::opType
  Ks::KsType
  KsCache::KsCacheType
end

function alg_cache(alg::ETDRK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if isa(f, SplitFunction)
    Jcache = nothing
  else
    # TODO: sparse Jacobian support
    Jcache = Matrix{eltype(u)}(length(u), length(u))
  end
  tmp = similar(u)
  rtmp = zeros(rate_prototype); Au = zeros(rate_prototype)
  F2 = zeros(rate_prototype); F3 = zeros(rate_prototype); F4 = zeros(rate_prototype)
  if alg.krylov
    ops = nothing # no caching
    n = length(u)
    m = min(alg.m, length(u))
    T = eltype(u)
    Ks = KrylovSubspace{T}(n, m)
    w1_half = Matrix{T}(n, 2); w2_half = Matrix{T}(n, 2)
    w1 = Matrix{T}(n, 4); w2 = Matrix{T}(n, 4); w3 = Matrix{T}(n, 4); w4 = Matrix{T}(n, 4)
    phiv_caches = (Vector{T}(m), Matrix{T}(m, m), Matrix{T}(m + 3, m + 3), Matrix{T}(m, 4))
    KsCache = (w1_half, w2_half, w1, w2, w3, w4, phiv_caches)
  else
    Ks = KsCache = nothing
    A = f.f1
    if isa(A, DiffEqArrayOperator)
      L = A.A .* A.α.coeff # has special handling if A.A is Diagonal
    else
      L = full(A)
    end
    ops = expRK_operators(alg, dt, L)
  end
  ETDRK4Cache(u,uprev,tmp,rtmp,Au,F2,F3,F4,Jcache,ops,Ks,KsCache)
end

u_cache(c::ETDRK4Cache) = ()
du_cache(c::ETDRK4Cache) = (c.k,c.fsalfirst,c.rtmp)

####################################
# Multistep exponential method caches

#=
  Fsal separately the linear and nonlinear part, as well as the nonlinear 
  part in the previous time step.
=#
mutable struct ETD2Fsal{rateType}
  lin::rateType
  nl::rateType
  nlprev::rateType
end
ETD2Fsal(rate_prototype) = ETD2Fsal(zero(rate_prototype),zero(rate_prototype),zero(rate_prototype))
function recursivecopy!(dest::ETD2Fsal, src::ETD2Fsal)
  recursivecopy!(dest.lin, src.lin)
  recursivecopy!(dest.nl, src.nl)
  recursivecopy!(dest.nlprev, src.nlprev)
end

struct ETD2ConstantCache{expType} <: OrdinaryDiffEqConstantCache
  exphA::expType
  phihA::expType
  B1::expType # ϕ1(hA) + ϕ2(hA)
  B0::expType # -ϕ2(hA)
end

function alg_cache(alg::ETD2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  A = f.f1
  if isa(A, DiffEqArrayOperator)
    _A = A.A * A.α.coeff # .* does not return Diagonal for A.A Diagonal
  else
    _A = full(A)
  end
  Phi = phi(dt*_A, 2)
  ETD2ConstantCache(Phi[1], Phi[2], Phi[2] + Phi[3], -Phi[3])
end

struct ETD2Cache{uType,rateType,expType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  utmp::uType
  rtmp1::rateType
  rtmp2::rateType
  exphA::expType
  phihA::expType
  B1::expType # ϕ1(hA) + ϕ2(hA)
  B0::expType # -ϕ2(hA)
end

function alg_cache(alg::ETD2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  A = f.f1
  if isa(A, DiffEqArrayOperator)
    _A = A.A * A.α.coeff # .* does not return Diagonal for A.A Diagonal
  else
    _A = full(A)
  end
  Phi = phi(dt*_A, 2)
  ETD2Cache(u,uprev,zero(u),zero(rate_prototype),zero(rate_prototype),Phi[1],Phi[2],Phi[2]+Phi[3],-Phi[3])
end

# TODO: what should these be?
u_cache(c::ETD2Cache) = ()
du_cache(c::ETD2Cache) = (c.rtmp1,c.rtmp2)
