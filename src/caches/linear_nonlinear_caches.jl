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
  uhold = Vector{typeof(u)}(undef, 1)
  rhs = RHS_IIF_Scalar(f,zero(u),t,t,one(uEltypeNoUnits),p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericIIF1ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::GenericIIF1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u,axes(u)); rtmp1 = zero(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  A = f.f1
  expA = exp(A*dt)
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
  uhold = Vector{typeof(u)}(undef, 1)
  tmp = zero(u)
  rhs = RHS_IIF_Scalar(f,tmp,t,t,uEltypeNoUnits(1//2),p)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  GenericIIF2ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::GenericIIF2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u,axes(u)); rtmp1 = zero(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  A = f.f1
  expA = exp(A*dt)
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
expRK_operators(::LawsonEuler, dt, A) = exp(dt * A)
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
function expRK_operators(::HochOst4, dt, A)
  P = phi(dt * A, 3)
  Phalf = phi(dt/2 * A, 3)
  A21 = 0.5Phalf[2]
  A31 = A21 - Phalf[3]
  A32 = Phalf[3]
  A41 = P[2] - 2P[3]
  A42 = P[3]
  A52 = 0.5Phalf[3] - P[4] + 0.25P[3] - 0.5Phalf[4]
  A54 = 0.25Phalf[3] - A52
  A51 = 0.5Phalf[2] - 2A52 - A54
  B1 = P[2] - 3P[3] + 4P[4]
  B4 = -P[3] + 4P[4]
  B5 = 4P[3] - 8P[4]
  return A21, A31, A32, A41, A42, A51, A52, A54, B1, B4, B5
end

# Unified constructor for constant caches
for (Alg, Cache) in [(:LawsonEuler, :LawsonEulerConstantCache),
                     (:NorsettEuler, :NorsettEulerConstantCache),
                     (:ETDRK2, :ETDRK2ConstantCache),
                     (:ETDRK3, :ETDRK3ConstantCache),
                     (:ETDRK4, :ETDRK4ConstantCache),
                     (:HochOst4, :HochOst4ConstantCache)]
  @eval struct $Cache{opType} <: ExpRKConstantCache
    ops::opType # precomputed operators
  end

  @eval function alg_cache(alg::$Alg,u,rate_prototype,uEltypeNoUnits,
    uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
    if alg.krylov
      ops = nothing # no caching
    else
      isa(f, SplitFunction) || throw(ArgumentError("Caching can only be used with SplitFunction"))
      A = size(f.f1) == () ? convert(Number, f.f1) : convert(AbstractMatrix, f.f1)
      ops = expRK_operators(alg, dt, A)
    end
    return $Cache(ops)
  end
end

"""
    alg_cache_expRK(alg,u,f,dt,plist)

Construct the non-standard caches (not uType or rateType) for ExpRK integrators.

`plist` is a list of integers each corresponding to the order of a `phiv(!)`
call in `perform_step!`.
"""
function alg_cache_expRK(alg::OrdinaryDiffEqExponentialAlgorithm, u, f, dt, plist)
  n = length(u); T = eltype(u)
  # Allocate cache for the Jacobian
  J = isa(f, SplitFunction) ? nothing : deepcopy(f.jac_prototype)
  if alg.krylov
    ops = nothing # no caching
    # Build up caches used by Krylov phiv
    m = min(alg.m, n)
    Ks = KrylovSubspace{T}(n, m)
    phiv_cache = PhivCache{T}(m, maximum(plist))
    ws = [Matrix{T}(undef, n, plist[i] + 1) for i = 1:length(plist)]
    KsCache = (Ks, phiv_cache, ws) # should use named tuple in v0.6
  else
    KsCache = nothing
    # Precompute the operators
    A = size(f.f1) == () ? convert(Number, f.f1) : convert(AbstractMatrix, f.f1)
    ops = expRK_operators(alg, dt, A)
  end
  return J, ops, KsCache
end

struct LawsonEulerCache{uType,rateType,JType,expType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  G::rateType
  J::JType
  exphA::expType
  KsCache::KsType
end

function alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                              # uType caches
  rtmp, G = (zero(rate_prototype) for i = 1:2)                 # rateType caches
  # other caches
  # This is different from other ExpRK integrators because LawsonEuler only
  # needs expv.
  n = length(u); T = eltype(u)
  J = isa(f, SplitFunction) ? nothing : deepcopy(f.jac_prototype)
  if alg.krylov
    exphA = nothing # no caching
    m = min(alg.m, n)
    Ks = KrylovSubspace{T}(n, m)
    expv_cache = ExpvCache{T}(m)
    KsCache = (Ks, expv_cache)
  else
    KsCache = nothing
    A = size(f.f1) == () ? convert(Number, f.f1) : convert(AbstractMatrix, f.f1)
    exphA = expRK_operators(alg, dt, A)
  end
  LawsonEulerCache(u,uprev,tmp,rtmp,G,J,exphA,KsCache)
end

u_cache(c::LawsonEulerCache) = ()
du_cache(c::LawsonEulerCache) = (c.rtmp)

struct NorsettEulerCache{uType,rateType,JType,expType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  G::rateType
  J::JType
  phihA::expType
  KsCache::KsType
end

function alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                              # uType caches
  rtmp, G = (zero(rate_prototype) for i = 1:2)                 # rateType caches
  J, phihA, KsCache = alg_cache_expRK(alg, u, f, dt, (1,)) # other caches
  NorsettEulerCache(u,uprev,tmp,rtmp,G,J,phihA,KsCache)
end

u_cache(c::NorsettEulerCache) = ()
du_cache(c::NorsettEulerCache) = (c.rtmp)

struct ETDRK2Cache{uType,rateType,JType,opType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  F2::rateType
  J::JType
  ops::opType
  KsCache::KsType
end

function alg_cache(alg::ETDRK2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                                  # uType caches
  rtmp, F2 = (zero(rate_prototype) for i = 1:2)                    # rateType caches
  J, ops, KsCache = alg_cache_expRK(alg, u, f, dt, (2, 2))   # other caches
  ETDRK2Cache(u,uprev,tmp,rtmp,F2,J,ops,KsCache)
end

struct ETDRK3Cache{uType,rateType,JType,opType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  Au::rateType
  F2::rateType
  F3::rateType
  J::JType
  ops::opType
  KsCache::KsType
end

function alg_cache(alg::ETDRK3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                                      # uType caches
  rtmp, Au, F2, F3 = (zero(rate_prototype) for i = 1:4)                # rateType caches
  J, ops, KsCache = alg_cache_expRK(alg, u, f, dt, (1,3,3,3))    # other caches
  ETDRK3Cache(u,uprev,tmp,rtmp,Au,F2,F3,J,ops,KsCache)
end

struct ETDRK4Cache{uType,rateType,JType,opType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  Au::rateType
  F2::rateType
  F3::rateType
  F4::rateType
  J::JType
  ops::opType
  KsCache::KsType
end

function alg_cache(alg::ETDRK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                                        # uType caches
  rtmp, Au, F2, F3, F4 = (zero(rate_prototype) for i = 1:5)              # rateType caches
  J, ops, KsCache = alg_cache_expRK(alg, u, f, dt, (1,1,3,3,3,3))  # other caches
  ETDRK4Cache(u,uprev,tmp,rtmp,Au,F2,F3,F4,J,ops,KsCache)
end

u_cache(c::ETDRK4Cache) = ()
du_cache(c::ETDRK4Cache) = (c.k,c.fsalfirst,c.rtmp)

struct HochOst4Cache{uType,rateType,JType,opType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  Au::rateType
  F2::rateType
  F3::rateType
  F4::rateType
  F5::rateType
  J::JType
  ops::opType
  KsCache::KsType
end

function alg_cache(alg::HochOst4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                                              # uType caches
  rtmp, rtmp2, Au, F2, F3, F4, F5 = (zero(rate_prototype) for i = 1:7)         # rateType caches
  J, ops, KsCache = alg_cache_expRK(alg, u, f, dt, (3,3,3,3,3,3,3,3,3))  # other caches
  HochOst4Cache(u,uprev,tmp,rtmp,rtmp2,Au,F2,F3,F4,F5,J,ops,KsCache)
end

####################################
# EPIRK method caches
struct Exp4ConstantCache <: ExpRKConstantCache end
alg_cache(alg::Exp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Exp4ConstantCache()

struct Exp4Cache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  K::matType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::Exp4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                              # uType caches
  rtmp, rtmp2 = (zero(rate_prototype) for i = 1:2)             # rateType caches
  # Allocate matrices
  # TODO: units
  n = length(u); T = eltype(u)
  K = Matrix{T}(undef, n, 3)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 2)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 1)
  Exp4Cache(u,uprev,tmp,rtmp,rtmp2,K,J,B,KsCache)
end

struct EPIRK4s3AConstantCache <: ExpRKConstantCache end
alg_cache(alg::EPIRK4s3A,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EPIRK4s3AConstantCache()

struct EPIRK4s3ACache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  K::matType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::EPIRK4s3A,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                    # uType caches
  rtmp, rtmp2 = (zero(rate_prototype) for i = 1:2)   # rateType caches
  # Allocate matrices
  n = length(u); T = eltype(u)
  K = Matrix{T}(undef, n, 2)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 5)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 4)
  EPIRK4s3ACache(u,uprev,tmp,rtmp,rtmp2,K,J,B,KsCache)
end

struct EPIRK4s3BConstantCache <: ExpRKConstantCache end
alg_cache(alg::EPIRK4s3B,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EPIRK4s3BConstantCache()

struct EPIRK4s3BCache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  K::matType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::EPIRK4s3B,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                    # uType caches
  rtmp, rtmp2 = (zero(rate_prototype) for i = 1:2)   # rateType caches
  # Allocate matrices
  n = length(u); T = eltype(u)
  K = Matrix{T}(undef, n, 2)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 5)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 4)
  EPIRK4s3BCache(u,uprev,tmp,rtmp,rtmp2,K,J,B,KsCache)
end

struct EPIRK5s3ConstantCache <: ExpRKConstantCache end
alg_cache(alg::EPIRK5s3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EPIRK5s3ConstantCache()

struct EPIRK5s3Cache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  k::uType
  rtmp::rateType
  rtmp2::rateType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::EPIRK5s3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp, k = (similar(u) for i = 1:2)                   # uType caches
  rtmp, rtmp2 = (zero(rate_prototype) for i = 1:2)   # rateType caches
  # Allocate matrices
  n = length(u); T = eltype(u)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 5)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 4)
  EPIRK5s3Cache(u,uprev,tmp,k,rtmp,rtmp2,J,B,KsCache)
end

struct EXPRB53s3ConstantCache <: ExpRKConstantCache end
alg_cache(alg::EXPRB53s3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EXPRB53s3ConstantCache()

struct EXPRB53s3Cache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  K::matType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::EXPRB53s3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                    # uType caches
  rtmp, rtmp2 = (zero(rate_prototype) for i = 1:2)   # rateType caches
  # Allocate matrices
  n = length(u); T = eltype(u)
  K = Matrix{T}(undef, n, 2)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 5)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 4)
  EXPRB53s3Cache(u,uprev,tmp,rtmp,rtmp2,K,J,B,KsCache)
end

struct EPIRK5P1ConstantCache <: ExpRKConstantCache end
alg_cache(alg::EPIRK5P1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EPIRK5P1ConstantCache()

struct EPIRK5P1Cache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  K::matType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::EPIRK5P1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                    # uType caches
  rtmp, rtmp2 = (zero(rate_prototype) for i = 1:2)   # rateType caches
  # Allocate matrices
  n = length(u); T = eltype(u)
  K = Matrix{T}(undef, n, 3)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 4)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 3)
  EPIRK5P1Cache(u,uprev,tmp,rtmp,rtmp2,K,J,B,KsCache)
end

struct EPIRK5P2ConstantCache <: ExpRKConstantCache end
alg_cache(alg::EPIRK5P2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
  uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EPIRK5P2ConstantCache()

struct EPIRK5P2Cache{uType,rateType,matType,JType,KsType} <: ExpRKCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
  rtmp2::rateType
  dR::rateType
  K::matType
  J::JType
  B::matType
  KsCache::KsType
end
function alg_cache(alg::EPIRK5P2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)                                        # uType caches
  rtmp, rtmp2, dR = (zero(rate_prototype) for i = 1:3)   # rateType caches
  # Allocate matrices
  n = length(u); T = eltype(u)
  K = Matrix{T}(undef, n, 3)
  J = deepcopy(f.jac_prototype)
  B = fill(zero(T), n, 4)
  # Allocate caches for phiv_timestep
  maxiter = min(alg.m, n)
  KsCache = _phiv_timestep_caches(u, maxiter, 3)
  EPIRK5P2Cache(u,uprev,tmp,rtmp,rtmp2,dR,K,J,B,KsCache)
end

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
  A = size(f.f1) == () ? convert(Number, f.f1) : convert(AbstractMatrix, f.f1)
  Phi = phi(dt*A, 2)
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
  A = size(f.f1) == () ? convert(Number, f.f1) : convert(AbstractMatrix, f.f1)
  Phi = phi(dt*A, 2)
  ETD2Cache(u,uprev,zero(u),zero(rate_prototype),zero(rate_prototype),Phi[1],Phi[2],Phi[2]+Phi[3],-Phi[3])
end

# TODO: what should these be?
u_cache(c::ETD2Cache) = ()
du_cache(c::ETD2Cache) = (c.rtmp1,c.rtmp2)
