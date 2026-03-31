################################################################################
# Implicit Weak Order 2 SRK Caches
# These caches are for drift-implicit versions of weak order 2 SRK methods

"""
    IRI1ConstantCache

Constant cache for the IRI1 (Implicit Rößler 1) method.
Contains the nlsolver for implicit drift treatment and the RI1 tableau coefficients.
"""
mutable struct IRI1ConstantCache{N, T, T2} <: StochasticDiffEqConstantCache
    nlsolver::N
    # RI1 Tableau coefficients (same as DRI1ConstantCache but stored directly)
    a021::T
    a031::T
    a032::T
    a121::T
    a131::T
    b021::T
    b031::T
    b121::T
    b131::T
    b221::T
    b222::T
    b223::T
    b231::T
    b232::T
    b233::T
    α1::T
    α2::T
    α3::T
    c02::T2
    c03::T2
    c12::T2
    c13::T2
    beta11::T
    beta12::T
    beta13::T
    beta22::T
    beta23::T
    beta31::T
    beta32::T
    beta33::T
    beta42::T
    beta43::T
    NORMAL_ONESIX_QUANTILE::T
end

function IRI1ConstantCache(nlsolver::N, ::Type{T}, ::Type{T2}) where {N, T, T2}
    # RI1 coefficients from the original RI1ConstantCache function
    a021 = convert(T, 2 // 3)
    a031 = convert(T, -1 // 3)
    a032 = convert(T, 1)

    a121 = convert(T, 1)
    a131 = convert(T, 1)

    b021 = convert(T, 1)
    b031 = convert(T, 0)

    b121 = convert(T, 1)
    b131 = convert(T, -1)

    b221 = convert(T, 1)
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -1)
    b232 = convert(T, 0)
    b233 = convert(T, 0)

    α1 = convert(T, 1 // 4)
    α2 = convert(T, 1 // 2)
    α3 = convert(T, 1 // 4)

    c02 = convert(T2, 2 // 3)
    c03 = convert(T2, 2 // 3)

    c12 = convert(T2, 1)
    c13 = convert(T2, 1)

    beta11 = convert(T, 1 // 2)
    beta12 = convert(T, 1 // 4)
    beta13 = convert(T, 1 // 4)

    beta22 = convert(T, 1 // 2)
    beta23 = convert(T, -1 // 2)

    beta31 = convert(T, -1 // 2)
    beta32 = convert(T, 1 // 4)
    beta33 = convert(T, 1 // 4)

    beta42 = convert(T, 1 // 2)
    beta43 = convert(T, -1 // 2)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return IRI1ConstantCache(
        nlsolver, a021, a031, a032, a121, a131, b021, b031, b121, b131,
        b221, b222, b223, b231, b232, b233, α1, α2, α3, c02, c03, c12, c13,
        beta11, beta12, beta13, beta22, beta23, beta31, beta32, beta33, beta42, beta43,
        NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::IRI1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = alg.theta, zero(t)
    nlsolver = OrdinaryDiffEqNonlinearSolve.build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return IRI1ConstantCache(nlsolver, real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

# IRI1Cache: Mutable cache for the IRI1 (Implicit Rößler 1) method (in-place version).
# Contains all the working arrays needed for the implicit weak order 2 SRK method.
@cache mutable struct IRI1Cache{
        uType, randType, rateNoiseType, rateType, noUnitsType, N, T, T2,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhat::uType
    _dW::randType
    _dZ::randType
    chi1::randType
    g1::rateNoiseType
    g2::Vector{rateNoiseType}
    g3::Vector{rateNoiseType}
    k1::rateType
    k2::rateType
    k3::rateType
    H02::uType
    H03::uType
    H12::Vector{uType}
    H13::Vector{uType}
    H22::Vector{uType}
    H23::Vector{uType}
    tmp::uType
    tmpg::rateNoiseType
    resids::noUnitsType
    nlsolver::N
    a021::T
    a031::T
    a032::T
    a121::T
    a131::T
    b021::T
    b031::T
    b121::T
    b131::T
    b221::T
    b222::T
    b223::T
    b231::T
    b232::T
    b233::T
    α1::T
    α2::T
    α3::T
    c02::T2
    c03::T2
    c12::T2
    c13::T2
    beta11::T
    beta12::T
    beta13::T
    beta22::T
    beta23::T
    beta31::T
    beta32::T
    beta33::T
    beta42::T
    beta43::T
    NORMAL_ONESIX_QUANTILE::T
end

function alg_cache(
        alg::IRI1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    T = real(uBottomEltypeNoUnits)
    T2 = real(tTypeNoUnits)

    # Build nonlinear solver
    γ, c = alg.theta, zero(t)
    nlsolver = OrdinaryDiffEqNonlinearSolve.build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )

    # Allocate arrays
    uhat = zero(u)

    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = ΔZ === nothing ? nothing : copy(ΔZ)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = ΔZ === nothing ? nothing : zero(ΔZ)
        chi1 = zero(ΔW)
    end

    m = is_diagonal_noise(prob) ? length(ΔW) : size(noise_rate_prototype, 2)

    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for _ in 1:m]
    g3 = [zero(noise_rate_prototype) for _ in 1:m]

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = [zero(u) for _ in 1:m]
    H13 = [zero(u) for _ in 1:m]
    H22 = [zero(u) for _ in 1:m]
    H23 = [zero(u) for _ in 1:m]

    tmp = zero(u)
    tmpg = zero(noise_rate_prototype)
    resids = fill!(similar(u, uEltypeNoUnits), zero(uEltypeNoUnits))

    # RI1 coefficients
    a021 = convert(T, 2 // 3)
    a031 = convert(T, -1 // 3)
    a032 = convert(T, 1)
    a121 = convert(T, 1)
    a131 = convert(T, 1)
    b021 = convert(T, 1)
    b031 = convert(T, 0)
    b121 = convert(T, 1)
    b131 = convert(T, -1)
    b221 = convert(T, 1)
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -1)
    b232 = convert(T, 0)
    b233 = convert(T, 0)
    α1 = convert(T, 1 // 4)
    α2 = convert(T, 1 // 2)
    α3 = convert(T, 1 // 4)
    c02 = convert(T2, 2 // 3)
    c03 = convert(T2, 2 // 3)
    c12 = convert(T2, 1)
    c13 = convert(T2, 1)
    beta11 = convert(T, 1 // 2)
    beta12 = convert(T, 1 // 4)
    beta13 = convert(T, 1 // 4)
    beta22 = convert(T, 1 // 2)
    beta23 = convert(T, -1 // 2)
    beta31 = convert(T, -1 // 2)
    beta32 = convert(T, 1 // 4)
    beta33 = convert(T, 1 // 4)
    beta42 = convert(T, 1 // 2)
    beta43 = convert(T, -1 // 2)
    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return IRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, g1, g2, g3, k1, k2, k3,
        H02, H03, H12, H13, H22, H23, tmp, tmpg, resids, nlsolver,
        a021, a031, a032, a121, a131, b021, b031, b121, b131,
        b221, b222, b223, b231, b232, b233, α1, α2, α3, c02, c03, c12, c13,
        beta11, beta12, beta13, beta22, beta23, beta31, beta32, beta33, beta42, beta43,
        NORMAL_ONESIX_QUANTILE
    )
end
