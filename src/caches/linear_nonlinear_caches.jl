#################################################
# Classical ExpRK method caches
abstract type ExpRKCache <: OrdinaryDiffEqMutableCache end
abstract type ExpRKConstantCache <: OrdinaryDiffEqConstantCache end

# Precomputation of exponential-like operators
"""
    expRK_operators(alg,dt,A) -> ops

Compute operator(s) for an ExpRK algorithm. `dt` is the time step and `A` is
the matrix form of the linear operator (from either a linear problem or a
SplitODEProblem). All ExpRK methods that use caching operators should implement
this method.
"""
function expRK_operators(alg::ExponentialAlgorithm, dt, A)
    error("$alg does not support caching operators at the moment.")
end
expRK_operators(::LawsonEuler, dt, A) = exp(dt * A)
expRK_operators(::NorsettEuler, dt, A) = phi(dt * A, 1)[2]
function expRK_operators(::ETDRK2, dt, A)
    P = phi(dt * A, 2)
    return P[2], P[3] # ϕ1(hA), ϕ2(hA)
end
function expRK_operators(::ETDRK3, dt, A)
    Phalf = phi(dt / 2 * A, 1)
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
    Phalf = phi(dt / 2 * A, 1)
    A21 = 0.5Phalf[2] # A32 = A21
    A41 = (dt / 4 * A) * Phalf[2]^2
    A43 = Phalf[2]
    B1 = P[2] - 3P[3] + 4P[4]
    B2 = 2P[3] - 4P[4] # B3 = B2
    B4 = -P[3] + 4P[4]
    return A21, A41, A43, B1, B2, B4
end
function expRK_operators(::HochOst4, dt, A)
    P = phi(dt * A, 3)
    Phalf = phi(dt / 2 * A, 3)
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
    @eval struct $Cache{opType, FType} <: ExpRKConstantCache
        ops::opType # precomputed operators
        uf::FType   # derivative wrapper
    end

    @eval function alg_cache(alg::$Alg, u, rate_prototype, ::Type{uEltypeNoUnits},
                             ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
                             uprev2, f, t, dt, reltol, p, calck,
                             ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits,
                                                  tTypeNoUnits}
        if alg.krylov
            ops = nothing # no caching
        else
            isa(f, SplitFunction) ||
                throw(ArgumentError("Caching can only be used with SplitFunction"))
            A = size(f.f1.f) == () ? convert(Number, f.f1.f) :
                convert(AbstractMatrix, f.f1.f)
            ops = expRK_operators(alg, dt, A)
        end
        if isa(f, SplitFunction) || DiffEqBase.has_jac(f)
            uf = nothing
        else
            uf = UDerivativeWrapper(f, t, p)
        end
        return $Cache(ops, uf)
    end
end

"""
    alg_cache_expRK(alg,u,uEltypeNoUnits,uprev,f,t,dt,p,du1,tmp,dz,plist)

Construct the non-standard caches (not uType or rateType) for ExpRK integrators.

`plist` is a list of integers each corresponding to the order of a `phiv(!)`
call in `perform_step!`.
"""
function alg_cache_expRK(alg::OrdinaryDiffEqExponentialAlgorithm, u, ::Type{uEltypeNoUnits},
                         uprev, f, t, dt, p, du1, tmp, dz, plist) where {uEltypeNoUnits}
    n = length(u)
    T = eltype(u)
    # Allocate cache for ForwardDiff
    if isa(f, SplitFunction) || DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    # Allocate cache for the Jacobian
    if isa(f, SplitFunction)
        J = nothing
    elseif DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), n, n)
    end
    if alg.krylov
        ops = nothing # no caching
        # Build up caches used by Krylov phiv
        m = min(alg.m, n)
        Ks = KrylovSubspace{T}(n, m)
        phiv_cache = PhivCache(u, m, maximum(plist))
        ws = [Matrix{T}(undef, n, plist[i] + 1) for i in 1:length(plist)]
        KsCache = (Ks, phiv_cache, ws)
    else
        KsCache = nothing
        # Precompute the operators
        A = size(f.f1.f) == () ? convert(Number, f.f1.f) : convert(AbstractMatrix, f.f1.f)
        ops = expRK_operators(alg, dt, A)
    end
    return uf, jac_config, J, ops, KsCache
end

@cache struct LawsonEulerCache{uType, rateType, JCType, FType, JType, expType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    G::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    exphA::expType
    KsCache::KsType
end

function alg_cache(alg::LawsonEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                            # uType caches
    rtmp, G, du1 = (zero(rate_prototype) for i in 1:3)             # rateType caches
    # other caches
    # This is different from other ExpRK integrators because LawsonEuler only
    # needs expv.
    n = length(u)
    T = eltype(u)
    # Allocate caches for ForwardDiff
    if isa(f, SplitFunction) || DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    # Allocate cache for the Jacobian
    if isa(f, SplitFunction)
        J = nothing
    elseif DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), n, n)
    end
    if alg.krylov
        exphA = nothing # no caching
        m = min(alg.m, n)
        Ks = KrylovSubspace{T}(n, m)
        expv_cache = ExpvCache{T}(m)
        KsCache = (Ks, expv_cache)
    else
        KsCache = nothing
        A = size(f.f1.f) == () ? convert(Number, f.f1.f) : convert(AbstractMatrix, f.f1.f)
        exphA = expRK_operators(alg, dt, A)
    end
    LawsonEulerCache(u, uprev, tmp, dz, rtmp, G, du1, jac_config, uf, J, exphA, KsCache)
end

@cache struct NorsettEulerCache{uType, rateType, JCType, FType, JType, expType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    G::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    phihA::expType
    KsCache::KsType
end

function alg_cache(alg::NorsettEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                # uType caches
    rtmp, G, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    plist = (1,)
    uf, jac_config, J, phihA, KsCache = alg_cache_expRK(alg, u, uEltypeNoUnits, uprev, f, t,
                                                        dt, p, du1, tmp, dz, plist) # other caches
    NorsettEulerCache(u, uprev, tmp, dz, rtmp, G, du1, jac_config, uf, J, phihA, KsCache)
end

@cache struct ETDRK2Cache{uType, rateType, JCType, FType, JType, opType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    F2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    ops::opType
    KsCache::KsType
end

function alg_cache(alg::ETDRK2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                 # uType caches
    rtmp, F2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    plist = (2, 2)
    uf, jac_config, J, ops, KsCache = alg_cache_expRK(alg, u, uEltypeNoUnits, uprev, f, t,
                                                      dt, p, du1, tmp, dz, plist) # other caches
    ETDRK2Cache(u, uprev, tmp, dz, rtmp, F2, du1, jac_config, uf, J, ops, KsCache)
end

@cache struct ETDRK3Cache{uType, rateType, JCType, FType, JType, opType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    Au::rateType
    F2::rateType
    F3::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    ops::opType
    KsCache::KsType
end

function alg_cache(alg::ETDRK3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                         # uType caches
    rtmp, Au, F2, F3, du1 = (zero(rate_prototype) for i in 1:5) # rateType caches
    plist = (1, 3, 3, 3)
    uf, jac_config, J, ops, KsCache = alg_cache_expRK(alg, u, uEltypeNoUnits, uprev, f, t,
                                                      dt, p, du1, tmp, dz, plist) # other caches
    ETDRK3Cache(u, uprev, tmp, dz, rtmp, Au, F2, F3, du1, jac_config, uf, J, ops, KsCache)
end

@cache struct ETDRK4Cache{uType, rateType, JCType, FType, JType, opType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    Au::rateType
    F2::rateType
    F3::rateType
    F4::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    ops::opType
    KsCache::KsType
end

function alg_cache(alg::ETDRK4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                             # uType caches
    rtmp, Au, F2, F3, F4, du1 = (zero(rate_prototype) for i in 1:6) # rateType caches
    plist = (1, 1, 3, 3, 3, 3)
    uf, jac_config, J, ops, KsCache = alg_cache_expRK(alg, u, uEltypeNoUnits, uprev, f, t,
                                                      dt, p, du1, tmp, dz, plist) # other caches
    ETDRK4Cache(u, uprev, tmp, dz, rtmp, Au, F2, F3, F4, du1, jac_config, uf, J, ops,
                KsCache)
end

@cache struct HochOst4Cache{uType, rateType, JCType, FType, JType, opType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    Au::rateType
    F2::rateType
    F3::rateType
    F4::rateType
    F5::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    ops::opType
    KsCache::KsType
end

function alg_cache(alg::HochOst4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                                        # uType caches
    rtmp, rtmp2, Au, F2, F3, F4, F5, du1 = (zero(rate_prototype) for i in 1:8) # rateType caches
    plist = (3, 3, 3, 3, 3, 3, 3, 3, 3)
    uf, jac_config, J, ops, KsCache = alg_cache_expRK(alg, u, uEltypeNoUnits, uprev, f, t,
                                                      dt, p, du1, tmp, dz, plist) # other caches
    HochOst4Cache(u, uprev, tmp, dz, rtmp, rtmp2, Au, F2, F3, F4, F5, du1, jac_config, uf,
                  J, ops, KsCache)
end

####################################
# EPIRK method caches
function _phiv_timestep_caches(u_prototype, maxiter::Int, p::Int)
    n = length(u_prototype)
    T = eltype(u_prototype)
    u = zero(u_prototype)                         # stores the current state
    W = Matrix{T}(undef, n, p + 1)                  # stores the w vectors
    P = Matrix{T}(undef, n, p + 2)                  # stores output from phiv!
    Ks = KrylovSubspace{T}(n, maxiter)            # stores output from arnoldi!
    phiv_cache = PhivCache(u_prototype, maxiter, p + 1) # cache used by phiv! (need +1 for error estimation)
    return u, W, P, Ks, phiv_cache
end

# Unified constructor for constant caches
for (Alg, Cache) in [(:Exp4, :Exp4ConstantCache),
    (:EPIRK4s3A, :EPIRK4s3AConstantCache),
    (:EPIRK4s3B, :EPIRK4s3BConstantCache),
    (:EPIRK5s3, :EPIRK5s3ConstantCache),
    (:EXPRB53s3, :EXPRB53s3ConstantCache),
    (:EPIRK5P1, :EPIRK5P1ConstantCache),
    (:EPIRK5P2, :EPIRK5P2ConstantCache)]
    @eval struct $Cache{FType} <: ExpRKConstantCache
        uf::FType   # derivative wrapper
    end
    @eval function alg_cache(alg::$Alg, u, rate_prototype, ::Type{uEltypeNoUnits},
                             ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
                             uprev2, f, t, dt, reltol, p, calck,
                             ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits,
                                                  tTypeNoUnits}
        if DiffEqBase.has_jac(f)
            uf = nothing
        else
            uf = UDerivativeWrapper(f, t, p)
        end
        return $Cache(uf)
    end
end

@cache struct Exp4Cache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    K::matType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::Exp4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                    # uType caches
    rtmp, rtmp2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    # TODO: units
    n = length(u)
    T = eltype(u)
    K = Matrix{T}(undef, n, 3)
    B = fill(zero(T), n, 2)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 1)
    Exp4Cache(u, uprev, tmp, dz, rtmp, rtmp2, du1, jac_config, uf, K, J, B, KsCache)
end

@cache struct EPIRK4s3ACache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    K::matType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::EPIRK4s3A, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                    # uType caches
    rtmp, rtmp2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    n = length(u)
    T = eltype(u)
    K = Matrix{T}(undef, n, 2)
    B = fill(zero(T), n, 5)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 4)
    EPIRK4s3ACache(u, uprev, tmp, dz, rtmp, rtmp2, du1, jac_config, uf, K, J, B, KsCache)
end

@cache struct EPIRK4s3BCache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    K::matType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::EPIRK4s3B, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                    # uType caches
    rtmp, rtmp2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    n = length(u)
    T = eltype(u)
    K = Matrix{T}(undef, n, 2)
    B = fill(zero(T), n, 5)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 4)
    EPIRK4s3BCache(u, uprev, tmp, dz, rtmp, rtmp2, du1, jac_config, uf, K, J, B, KsCache)
end

@cache struct EPIRK5s3Cache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    k::uType
    rtmp::rateType
    rtmp2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::EPIRK5s3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz, k = (zero(u) for i in 1:3)                 # uType caches
    rtmp, rtmp2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    n = length(u)
    T = eltype(u)
    B = fill(zero(T), n, 5)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 4)
    EPIRK5s3Cache(u, uprev, tmp, dz, k, rtmp, rtmp2, du1, jac_config, uf, J, B, KsCache)
end

@cache struct EXPRB53s3Cache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    K::matType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::EXPRB53s3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                    # uType caches
    rtmp, rtmp2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    n = length(u)
    T = eltype(u)
    K = Matrix{T}(undef, n, 2)
    B = fill(zero(T), n, 5)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 4)
    EXPRB53s3Cache(u, uprev, tmp, dz, rtmp, rtmp2, du1, jac_config, uf, K, J, B, KsCache)
end

@cache struct EPIRK5P1Cache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    K::matType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::EPIRK5P1, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                    # uType caches
    rtmp, rtmp2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    n = length(u)
    T = eltype(u)
    K = Matrix{T}(undef, n, 3)
    B = fill(zero(T), n, 4)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 3)
    EPIRK5P1Cache(u, uprev, tmp, dz, rtmp, rtmp2, du1, jac_config, uf, K, J, B, KsCache)
end

@cache struct EPIRK5P2Cache{uType, rateType, JCType, FType, matType, JType, KsType} <:
              ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    rtmp2::rateType
    dR::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    K::matType
    J::JType
    B::matType
    KsCache::KsType
end
function alg_cache(alg::EPIRK5P2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits},
                   ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp, dz = (zero(u) for i in 1:2)                        # uType caches
    rtmp, rtmp2, dR, du1 = (zero(rate_prototype) for i in 1:4) # rateType caches
    # Allocate jacobian and caches for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), length(u), length(u))
    end
    # Allocate matrices
    n = length(u)
    T = eltype(u)
    K = Matrix{T}(undef, n, 3)
    B = fill(zero(T), n, 4)
    # Allocate caches for phiv_timestep
    maxiter = min(alg.m, n)
    KsCache = _phiv_timestep_caches(u, maxiter, 3)
    EPIRK5P2Cache(u, uprev, tmp, dz, rtmp, rtmp2, dR, du1, jac_config, uf, K, J, B, KsCache)
end

####################################
# Adaptive exponential Rosenbrock method caches

## Constant caches
for (Alg, Cache) in [(:Exprb32, :Exprb32ConstantCache),
    (:Exprb43, :Exprb43ConstantCache)]
    @eval struct $Cache{FType} <: ExpRKConstantCache
        uf::FType   # derivative wrapper
    end
    @eval function alg_cache(alg::$Alg, u, rate_prototype, ::Type{uEltypeNoUnits},
                             ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
                             uprev2, f, t, dt, reltol, p, calck,
                             ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits,
                                                  tTypeNoUnits}
        if DiffEqBase.has_jac(f)
            uf = nothing
        else
            uf = UDerivativeWrapper(f, t, p)
        end
        return $Cache(uf)
    end
end

## Mutable caches
"""
    alg_cache_exprb(alg,uEltypeNoUnits,uprev,f,t,p,du1,tmp,dz,plist)

Construct the non-standard caches (not uType or rateType) for Exprb integrators.

`plist` is a list of integers each corresponding to the order of a `phiv(!)`
call in `perform_step!`.
"""
function alg_cache_exprb(alg::OrdinaryDiffEqAdaptiveExponentialAlgorithm, u,
                         ::Type{uEltypeNoUnits}, uprev, f, t, p, du1, tmp, dz,
                         plist) where {uEltypeNoUnits}
    if f isa SplitFunction
        error("Algorithm $alg cannot be used for split problems. Consider reformat to a regular `ODEProblem`")
    end
    n = length(u)
    T = eltype(u)
    # Allocate cache for ForwardDiff
    if DiffEqBase.has_jac(f)
        uf = nothing
        jac_config = nothing
    else
        uf = UJacobianWrapper(f, t, p)
        jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dz)
    end
    # Allocate cache for the Jacobian
    if DiffEqBase.has_jac(f) && f.jac_prototype !== nothing
        J = deepcopy(f.jac_prototype)
    else
        J = fill(zero(uEltypeNoUnits), n, n)
    end
    # Build up caches used by Krylov phiv
    m = min(alg.m, n)
    Ks = KrylovSubspace{T}(n, m)
    phiv_cache = PhivCache(u, m, maximum(plist))
    ws = [Matrix{T}(undef, n, plist[i] + 1) for i in 1:length(plist)]
    KsCache = (Ks, phiv_cache, ws)
    return uf, jac_config, J, KsCache
end

@cache struct Exprb32Cache{uType, rateType, JCType, FType, JType, KsType} <: ExpRKCache
    u::uType
    uprev::uType
    utilde::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    F2::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    KsCache::KsType
end
function alg_cache(alg::Exprb32, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utilde, tmp, dz = (zero(u) for i in 1:3)         # uType caches
    rtmp, F2, du1 = (zero(rate_prototype) for i in 1:3) # rateType caches
    plist = (3, 3)
    uf, jac_config, J, KsCache = alg_cache_exprb(alg, u, uEltypeNoUnits, uprev, f, t, p,
                                                 du1, tmp, dz, plist) # other caches
    Exprb32Cache(u, uprev, utilde, tmp, dz, rtmp, F2, du1, jac_config, uf, J, KsCache)
end

struct Exprb43Cache{uType, rateType, JCType, FType, JType, KsType} <: ExpRKCache
    u::uType
    uprev::uType
    utilde::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    Au::rateType
    F2::rateType
    F3::rateType
    du1::rateType
    jac_config::JCType
    uf::FType
    J::JType
    KsCache::KsType
end
function alg_cache(alg::Exprb43, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utilde, tmp, dz = (zero(u) for i in 1:3)                 # uType caches
    rtmp, Au, F2, F3, du1 = (zero(rate_prototype) for i in 1:5) # rateType caches
    plist = (1, 4, 4, 4)
    uf, jac_config, J, KsCache = alg_cache_exprb(alg, u, uEltypeNoUnits, uprev, f, t, p,
                                                 du1, tmp, dz, plist) # other caches
    Exprb43Cache(u, uprev, utilde, tmp, dz, rtmp, Au, F2, F3, du1, jac_config, uf, J,
                 KsCache)
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

    function ETD2Fsal(lin, nl, nlprev)
        if size(lin) == ()
            # convert to same type if Number or AbstractSciMLScalarOperator
            T = promote_type(eltype.((lin, nl, nlprev))...)

            lin = convert(T, lin)
            nl = convert(T, nl)
            nlprev = convert(T, nlprev)
        end

        new{typeof(lin)}(lin, nl, nlprev)
    end
end
function ETD2Fsal(rate_prototype)
    ETD2Fsal(zero(rate_prototype), zero(rate_prototype), zero(rate_prototype))
end
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

function alg_cache(alg::ETD2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    A = size(f.f1.f) == () ? convert(Number, f.f1.f) : convert(AbstractMatrix, f.f1.f)
    Phi = phi(dt * A, 2)
    ETD2ConstantCache(Phi[1], Phi[2], Phi[2] + Phi[3], -Phi[3])
end

@cache struct ETD2Cache{uType, rateType, expType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::ETD2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    A = size(f.f1.f) == () ? convert(Number, f.f1.f) : convert(AbstractMatrix, f.f1.f)
    Phi = phi(dt * A, 2)
    ETD2Cache(u, uprev, zero(u), zero(rate_prototype), zero(rate_prototype), Phi[1], Phi[2],
              Phi[2] + Phi[3], -Phi[3])
end
