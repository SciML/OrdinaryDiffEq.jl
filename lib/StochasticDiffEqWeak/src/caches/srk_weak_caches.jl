################################################################################
# Roessler SRK for second order weak approx

struct DRI1ConstantCache{T, T2} <: StochasticDiffEqConstantCache
    # hard-coded version
    a021::T
    a031::T
    a032::T

    a121::T
    a131::T
    #a132::T

    #a221::T
    #a231::T
    #a232::T

    b021::T
    b031::T
    #b032::T

    b121::T
    b131::T
    #b132::T

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

    #c11::T2
    c12::T2
    c13::T2

    #c21::T2
    #c22::T2
    #c23::T2

    beta11::T
    beta12::T
    beta13::T

    #beta21::T
    beta22::T
    beta23::T

    beta31::T
    beta32::T
    beta33::T

    #beta41::T
    beta42::T
    beta43::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function DRI1ConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    a021 = convert(T, 1 // 2)
    a031 = convert(T, -1)
    a032 = convert(T, 2)

    a121 = convert(T, 342 // 491)
    a131 = convert(T, 342 // 491)

    b021 = convert(T, (6 - sqrt(6)) / 10)
    b031 = convert(T, (3 + 2 * sqrt(6)) / 5)

    b121 = convert(T, 3 * sqrt(38 // 491))
    b131 = convert(T, -3 * sqrt(38 / 491))

    b221 = convert(T, -214 // 513 * sqrt(1105 // 991))
    b222 = convert(T, -491 // 513 * sqrt(221 // 4955))
    b223 = convert(T, -491 // 513 * sqrt(221 // 4955))
    b231 = convert(T, 214 // 513 * sqrt(1105 // 991))
    b232 = convert(T, 491 // 513 * sqrt(221 // 4955))
    b233 = convert(T, 491 // 513 * sqrt(221 // 4955))

    α1 = convert(T, 1 // 6)
    α2 = convert(T, 2 // 3)
    α3 = convert(T, 1 // 6)

    c02 = convert(T2, 1 // 2)
    c03 = convert(T2, 1)

    c12 = convert(T2, 342 // 491)
    c13 = convert(T2, 342 // 491)

    beta11 = convert(T, 193 // 684)
    beta12 = convert(T, 491 // 1368)
    beta13 = convert(T, 491 // 1368)

    beta22 = convert(T, 1 // 6 * sqrt(491 // 38))
    beta23 = convert(T, -1 // 6 * sqrt(491 // 38))

    beta31 = convert(T, -4955 // 7072)
    beta32 = convert(T, 4955 // 14144)
    beta33 = convert(T, 4955 // 14144)

    beta42 = convert(T, -1 // 8 * sqrt(4955 // 221))
    beta43 = convert(T, 1 // 8 * sqrt(4955 // 221))

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::DRI1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RI1ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 2 // 3) # convert(T, 2//3)
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

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RI1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RI3ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1)
    a031 = convert(T, 1 // 4)
    a032 = convert(T, 1 // 4)

    a121 = convert(T, 1)
    a131 = convert(T, 1)

    b021 = convert(T, (3 - 2 * sqrt(6)) / 5)
    b031 = convert(T, (6 + sqrt(6)) / 10)

    b121 = convert(T, 1)
    b131 = convert(T, -1)

    b221 = convert(T, 1)
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -1)
    b232 = convert(T, 0)
    b233 = convert(T, 0)

    α1 = convert(T, 1 // 6)
    α2 = convert(T, 1 // 6)
    α3 = convert(T, 2 // 3)

    c02 = convert(T2, 1)
    c03 = convert(T2, 1 // 2)

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

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RI3, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RI3ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RI5ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1)
    a031 = convert(T, 25 // 144)
    a032 = convert(T, 35 // 144)

    a121 = convert(T, 1 // 4)
    a131 = convert(T, 1 // 4)

    b021 = convert(T, 1 // 3)
    b031 = convert(T, -5 // 6)

    b121 = convert(T, 1 // 2)
    b131 = convert(T, -1 // 2)

    b221 = convert(T, 1)
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -1)
    b232 = convert(T, 0)
    b233 = convert(T, 0)

    α1 = convert(T, 1 // 10)
    α2 = convert(T, 3 // 14)
    α3 = convert(T, 24 // 35)

    c02 = convert(T2, 1)
    c03 = convert(T2, 5 // 12)

    c12 = convert(T2, 1 // 4)
    c13 = convert(T2, 1 // 4)

    beta11 = convert(T, 1)
    beta12 = convert(T, -1)
    beta13 = convert(T, -1)

    beta22 = convert(T, 1)
    beta23 = convert(T, -1)

    beta31 = convert(T, 1 // 2)
    beta32 = convert(T, -1 // 4)
    beta33 = convert(T, -1 // 4)

    beta42 = convert(T, 1 // 2)
    beta43 = convert(T, -1 // 2)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RI5, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RI5ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RI6ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1)
    a031 = convert(T, 0)
    a032 = convert(T, 0)

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

    α1 = convert(T, 1 // 2)
    α2 = convert(T, 1 // 2)
    α3 = convert(T, 0)

    c02 = convert(T2, 1)
    c03 = convert(T2, 0)

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

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RI6, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RI6ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RDI2WMConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1)
    a031 = convert(T, 0)
    a032 = convert(T, 0)

    a121 = convert(T, 2 // 3)
    a131 = convert(T, 2 // 3)

    b021 = convert(T, 1)
    b031 = convert(T, 0)

    b121 = convert(T, sqrt(2 // 3))
    b131 = convert(T, -sqrt(2 // 3))

    b221 = convert(T, sqrt(2))
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -sqrt(2))
    b232 = convert(T, 0)
    b233 = convert(T, 0)

    α1 = convert(T, 1 // 2)
    α2 = convert(T, 1 // 2)
    α3 = convert(T, 0)

    c02 = convert(T2, 1)
    c03 = convert(T2, 0)

    c12 = convert(T2, 2 // 3)
    c13 = convert(T2, 2 // 3)

    beta11 = convert(T, 1 // 4)
    beta12 = convert(T, 3 // 8)
    beta13 = convert(T, 3 // 8)

    beta22 = convert(T, sqrt(6) / 4)
    beta23 = convert(T, -sqrt(6) / 4)

    beta31 = convert(T, -1 // 4)
    beta32 = convert(T, 1 // 8)
    beta33 = convert(T, 1 // 8)

    beta42 = convert(T, sqrt(2) / 4)
    beta43 = convert(T, -sqrt(2) / 4)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RDI2WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RDI2WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RDI3WMConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1 // 2)
    a031 = convert(T, 0)
    a032 = convert(T, 3 // 4)

    a121 = convert(T, 2 // 3)
    a131 = convert(T, 2 // 3)

    b021 = convert(T, (9 - 2 * sqrt(15)) / 14)
    b031 = convert(T, (18 + 3 * sqrt(15)) / 28)

    b121 = convert(T, sqrt(2 // 3))
    b131 = convert(T, -sqrt(2 // 3))

    b221 = convert(T, sqrt(2))
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -sqrt(2))
    b232 = convert(T, 0)
    b233 = convert(T, 0)

    α1 = convert(T, 2 // 9)
    α2 = convert(T, 1 // 3)
    α3 = convert(T, 4 // 9)

    c02 = convert(T2, 1 // 2)
    c03 = convert(T2, 3 // 4)

    c12 = convert(T2, 2 // 3)
    c13 = convert(T2, 2 // 3)

    beta11 = convert(T, 1 // 4)
    beta12 = convert(T, 3 // 8)
    beta13 = convert(T, 3 // 8)

    beta22 = convert(T, sqrt(6) / 4)
    beta23 = convert(T, -sqrt(6) / 4)

    beta31 = convert(T, -1 // 4)
    beta32 = convert(T, 1 // 8)
    beta33 = convert(T, 1 // 8)

    beta42 = convert(T, sqrt(2) / 4)
    beta43 = convert(T, -sqrt(2) / 4)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RDI3WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RDI3WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RDI4WMConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1 // 2)
    a031 = convert(T, -1)
    a032 = convert(T, 2)

    a121 = convert(T, 2 // 3)
    a131 = convert(T, 2 // 3)

    b021 = convert(T, (6 - sqrt(6)) / 10)
    b031 = convert(T, (3 + 2 * sqrt(6)) / 5)

    b121 = convert(T, sqrt(2 // 3))
    b131 = convert(T, -sqrt(2 // 3))

    b221 = convert(T, sqrt(2))
    b222 = convert(T, 0)
    b223 = convert(T, 0)
    b231 = convert(T, -sqrt(2))
    b232 = convert(T, 0)
    b233 = convert(T, 0)

    α1 = convert(T, 1 // 6)
    α2 = convert(T, 2 // 3)
    α3 = convert(T, 1 // 6)

    c02 = convert(T2, 1 // 2)
    c03 = convert(T2, 1)

    c12 = convert(T2, 2 // 3)
    c13 = convert(T2, 2 // 3)

    beta11 = convert(T, 1 // 4)
    beta12 = convert(T, 3 // 8)
    beta13 = convert(T, 3 // 8)

    beta22 = convert(T, sqrt(6) / 4)
    beta23 = convert(T, -sqrt(6) / 4)

    beta31 = convert(T, -1 // 4)
    beta32 = convert(T, 1 // 8)
    beta33 = convert(T, 1 // 8)

    beta42 = convert(T, sqrt(2) / 4)
    beta43 = convert(T, -sqrt(2) / 4)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return DRI1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, b222, b223, b231,
        b232, b233, α1, α2, α3, c02, c03, c12, c13, beta11, beta12, beta13, beta22,
        beta23, beta31, beta32, beta33, beta42, beta43, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RDI4WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RDI4WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct DRI1Cache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhat::uType

    _dW::randType
    _dZ::randType
    chi1::randType

    tab::tabType

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

    tmp1::possibleRateType
    tmpg::rateNoiseType

    tmp::uType
    resids::uType
end

@cache struct DRI1NMCache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhat::uType

    _dW::randType
    chi1::randType

    tab::tabType

    g1::rateNoiseType
    g2::rateNoiseType
    g3::rateNoiseType

    k1::rateType
    k2::rateType
    k3::rateType

    H02::uType
    H03::uType
    H12::rateNoiseType
    H13::rateNoiseType

    tmp1::possibleRateType
    tmpg::rateNoiseType

    tmp::uType
    resids::uType
end

function alg_cache(
        alg::DRI1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::DRI1NM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = zero(noise_rate_prototype)
    H13 = zero(noise_rate_prototype)
    #H22 = zero(noise_rate_prototype)
    #H23 = zero(noise_rate_prototype)

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1NMCache(
        u, uprev, uhat, _dW, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RI1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RI3, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RI3ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RI5, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RI5ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RI6, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RI6ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RDI2WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RDI2WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RDI3WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RDI3WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

function alg_cache(
        alg::RDI4WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RDI4WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return DRI1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1, k2,
        k3, H02, H03, H12, H13, H22, H23, tmp1, tmpg, tmp, resids
    )
end

# Roessler SRK for first order weak approx
struct RDI1WMConstantCache{T, T2} <: StochasticDiffEqConstantCache
    # hard-coded version
    a021::T

    #a121::T

    #a221::T

    b021::T

    #b121::T

    #b221::T

    α1::T
    α2::T

    c02::T2

    #c11::T2
    #c12::T2

    #c21::T2
    #c22::T2

    beta11::T
    #beta12::T

    #beta21::T
    #beta22::T

    #beta31::T
    #beta32::T

    #beta41::T
    #beta42::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function RDI1WMConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    a021 = convert(T, 2 // 3)

    b021 = convert(T, 2 // 3)

    α1 = convert(T, 1 // 4)
    α2 = convert(T, 3 // 4)

    c02 = convert(T2, 2 // 3)

    beta11 = convert(T, 1)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return RDI1WMConstantCache(a021, b021, α1, α2, c02, beta11, NORMAL_ONESIX_QUANTILE)
end

function alg_cache(
        alg::RDI1WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RDI1WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct RDI1WMCache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType
    chi1::randType

    tab::tabType

    g1::rateNoiseType

    k1::rateType
    k2::rateType

    H02::uType

    tmp1::possibleRateType
end

function alg_cache(
        alg::RDI1WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RDI1WMConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)

    H02 = zero(u)

    tmp1 = zero(rate_prototype)

    return RDI1WMCache(u, uprev, _dW, chi1, tab, g1, k1, k2, H02, tmp1)
end

# Stratonovich sense

struct RSConstantCache{T, T2} <: StochasticDiffEqConstantCache
    # hard-coded version
    a021::T
    a031::T
    a032::T

    a131::T
    a141::T

    b031::T
    b032::T

    b121::T
    b131::T
    b132::T
    b141::T
    b142::T
    b143::T

    b221::T
    b231::T

    b331::T
    b332::T
    b341::T
    b342::T

    α1::T
    α2::T
    α3::T
    α4::T

    c02::T2
    c03::T2

    c13::T2
    c14::T2

    beta11::T
    beta12::T
    beta13::T
    beta14::T

    beta22::T
    beta23::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function RS1ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 0)
    a031 = convert(T, 1)
    a032 = convert(T, 0)

    a131 = convert(T, 1)
    a141 = convert(T, 1)

    b031 = convert(T, 1 // 4)
    b032 = convert(T, 3 // 4)

    b121 = convert(T, 2 // 3)
    b131 = convert(T, 1 // 12)
    b132 = convert(T, 1 // 4)
    b141 = convert(T, -5 // 4)
    b142 = convert(T, 1 // 4)
    b143 = convert(T, 2)

    b221 = convert(T, 1)
    b231 = convert(T, -1)

    b331 = convert(T, 1 // 4)
    b332 = convert(T, 3 // 4)
    b341 = convert(T, 1 // 4)
    b342 = convert(T, 3 // 4)

    α1 = convert(T, 0)
    α2 = convert(T, 0)
    α3 = convert(T, 1 // 2)
    α4 = convert(T, 1 // 2)

    c02 = convert(T2, 0)
    c03 = convert(T2, 1)

    c13 = convert(T2, 1)
    c14 = convert(T2, 1)

    beta11 = convert(T, 1 // 8)
    beta12 = convert(T, 3 // 8)
    beta13 = convert(T, 3 // 8)
    beta14 = convert(T, 1 // 8)

    beta22 = convert(T, -1 // 4)
    beta23 = convert(T, 1 // 4)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return RSConstantCache(
        a021, a031, a032, a131, a141, b031, b032, b121, b131, b132, b141, b142, b143,
        b221, b231, b331, b332, b341, b342, α1, α2, α3, α4, c02, c03, c13, c14,
        beta11, beta12, beta13, beta14, beta22, beta23, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RS1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RS1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function RS2ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 2 // 3)
    a031 = convert(T, 1 // 6)
    a032 = convert(T, 1 // 2)

    a131 = convert(T, 1)
    a141 = convert(T, 1)

    b031 = convert(T, 1 // 4)
    b032 = convert(T, 3 // 4)

    b121 = convert(T, 2 // 3)
    b131 = convert(T, 1 // 12)
    b132 = convert(T, 1 // 4)
    b141 = convert(T, -5 // 4)
    b142 = convert(T, 1 // 4)
    b143 = convert(T, 2)

    b221 = convert(T, 1)
    b231 = convert(T, -1)

    b331 = convert(T, 1 // 4)
    b332 = convert(T, 3 // 4)
    b341 = convert(T, 1 // 4)
    b342 = convert(T, 3 // 4)

    α1 = convert(T, 1 // 4)
    α2 = convert(T, 1 // 4)
    α3 = convert(T, 1 // 2)
    α4 = convert(T, 0)

    c02 = convert(T2, 2 // 3)
    c03 = convert(T2, 2 // 3)

    c13 = convert(T2, 1)
    c14 = convert(T2, 1)

    beta11 = convert(T, 1 // 8)
    beta12 = convert(T, 3 // 8)
    beta13 = convert(T, 3 // 8)
    beta14 = convert(T, 1 // 8)

    beta22 = convert(T, -1 // 4)
    beta23 = convert(T, 1 // 4)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return RSConstantCache(
        a021, a031, a032, a131, a141, b031, b032, b121, b131, b132, b141, b142, b143,
        b221, b231, b331, b332, b341, b342, α1, α2, α3, α4, c02, c03, c13, c14,
        beta11, beta12, beta13, beta14, beta22, beta23, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::RS2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RS2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct RSCache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType
    _dZ::randType
    chi1::randType

    tab::tabType

    g1::rateNoiseType
    g2::Vector{rateNoiseType}
    g3::Vector{rateNoiseType}
    g4::Vector{rateNoiseType}

    k1::rateType
    k2::rateType
    k3::rateType

    H02::uType
    H03::uType
    H12::Vector{uType}
    H13::Vector{uType}
    H14::Vector{uType}
    H22::Vector{uType}
    H23::Vector{uType}

    tmp1::possibleRateType
    tmpg::rateNoiseType
end

function alg_cache(
        alg::RS1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RS1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    g4 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H14 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H14, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    return RSCache(
        u, uprev, _dW, _dZ, chi1, tab, g1, g2, g3, g4, k1, k2,
        k3, H02, H03, H12, H13, H14, H22, H23, tmp1, tmpg
    )
end

function alg_cache(
        alg::RS2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = RS2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = [zero(noise_rate_prototype) for k in 1:m]
    g3 = [zero(noise_rate_prototype) for k in 1:m]
    g4 = [zero(noise_rate_prototype) for k in 1:m]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()
    H14 = Vector{typeof(u)}()
    H22 = Vector{typeof(u)}()
    H23 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
        push!(H14, zero(u))
        push!(H22, zero(u))
        push!(H23, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    return RSCache(
        u, uprev, _dW, _dZ, chi1, tab, g1, g2, g3, g4, k1, k2,
        k3, H02, H03, H12, H13, H14, H22, H23, tmp1, tmpg
    )
end

# PL1WM
struct PL1WMConstantCache{T} <: StochasticDiffEqConstantCache
    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function PL1WMConstantCache(T::Type)
    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return PL1WMConstantCache(NORMAL_ONESIX_QUANTILE)
end

function alg_cache(
        alg::PL1WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return PL1WMConstantCache(real(uBottomEltypeNoUnits))
end

struct PL1WMAConstantCache{T} <: StochasticDiffEqConstantCache
    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function PL1WMAConstantCache(T::Type)
    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return PL1WMAConstantCache(NORMAL_ONESIX_QUANTILE)
end

function alg_cache(
        alg::PL1WMA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return PL1WMAConstantCache(real(uBottomEltypeNoUnits))
end

@cache struct PL1WMCache{
        uType, randType, rand2Type, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType
    _dZ::rand2Type
    chi1::randType

    tab::tabType

    g1::rateNoiseType

    k1::rateType
    k2::rateType

    Y::uType
    Yp::Vector{uType}
    Ym::Vector{uType}

    tmp1::possibleRateType
    tmpg1::rateNoiseType
    tmpg2::rateNoiseType
    Ulp::uType
    Ulm::uType
end

function alg_cache(
        alg::PL1WM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔZ)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔZ)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = PL1WMConstantCache(real(uBottomEltypeNoUnits))
    g1 = zero(noise_rate_prototype)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)

    Y = zero(u)
    Yp = Vector{typeof(u)}()
    Ym = Vector{typeof(u)}()

    for k in 1:m
        push!(Yp, zero(u))
        push!(Ym, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg1 = zero(noise_rate_prototype)
    tmpg2 = zero(noise_rate_prototype)

    Ulp = zero(u)
    Ulm = zero(u)

    return PL1WMCache(
        u, uprev, _dW, _dZ, chi1, tab, g1, k1, k2, Y, Yp, Ym, tmp1, tmpg1, tmpg2, Ulp, Ulm
    )
end

# additive noise
@cache struct PL1WMACache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType
    chi1::randType
    tab::tabType

    g1::rateNoiseType
    k1::rateType
    k2::rateType

    Y::uType
    tmp1::possibleRateType
end

function alg_cache(
        alg::PL1WMA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        chi1 = zero(ΔW)
    end

    tab = PL1WMConstantCache(real(uBottomEltypeNoUnits))
    g1 = zero(noise_rate_prototype)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)

    Y = zero(u)

    tmp1 = zero(rate_prototype)

    return PL1WMACache(u, uprev, _dW, chi1, tab, g1, k1, k2, Y, tmp1)
end

#NON
struct NONConstantCache{T} <: StochasticDiffEqConstantCache
    c01::T
    c02::T
    c03::T
    c04::T

    cj1::T
    cj2::T
    cj3::T
    cj4::T

    cjl2::T
    cjl3::T

    clj2::T
    clj3::T

    a0021::T
    a0032::T
    a0043::T

    aj021::T
    aj041::T

    a0j21::T
    a0j31::T
    a0j32::T
    a0j41::T

    ajj21::T
    ajj31::T
    ajj32::T
    ajj41::T
    ajj42::T
    ajj43::T

    ajl31::T
    ajl32::T
    ajl41::T
    ajl42::T

    ajljj31::T

    aljjl21::T
    aljjl31::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function NONConstantCache(T::Type)
    c01 = convert(T, 1 // 6)
    c02 = convert(T, 1 // 3)
    c03 = convert(T, 1 // 3)
    c04 = convert(T, 1 // 6)

    cj1 = convert(T, 1 // 8)
    cj2 = convert(T, 3 // 8)
    cj3 = convert(T, 3 // 8)
    cj4 = convert(T, 1 // 8)

    cjl2 = convert(T, 1 // 2)
    cjl3 = convert(T, -1 // 2)

    clj2 = convert(T, -1 // 2)
    clj3 = convert(T, 1 // 2)

    a0021 = convert(T, 1 // 2)
    a0032 = convert(T, 1 // 2)
    a0043 = convert(T, 1)

    aj021 = convert(T, 2)
    aj041 = convert(T, -2)

    a0j21 = convert(T, 1)
    a0j31 = convert(T, -9 // 8)
    a0j32 = convert(T, 9 // 8)
    a0j41 = convert(T, 1)

    ajj21 = convert(T, 2 // 3)
    ajj31 = convert(T, 1 // 12)
    ajj32 = convert(T, 1 // 4)
    ajj41 = convert(T, -5 // 4)
    ajj42 = convert(T, 1 // 4)
    ajj43 = convert(T, 2)

    ajl31 = convert(T, 1 // 4)
    ajl32 = convert(T, 3 // 4)
    ajl41 = convert(T, 1 // 4)
    ajl42 = convert(T, 3 // 4)

    ajljj31 = convert(T, 1)

    aljjl21 = convert(T, -1 // 2)
    aljjl31 = convert(T, 1 // 2)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return NONConstantCache(
        c01, c02, c03, c04, cj1, cj2, cj3, cj4, cjl2, cjl3, clj2, clj3, a0021, a0032, a0043,
        aj021, aj041, a0j21, a0j31, a0j32, a0j41, ajj21, ajj31, ajj32, ajj41, ajj42, ajj43,
        ajl31, ajl32, ajl41, ajl42, ajljj31, aljjl21, aljjl31, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::NON, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NONConstantCache(real(uBottomEltypeNoUnits))
end

@cache struct NONCache{uType, randType, tabType, rateNoiseType, rateType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType
    _dZ::randType
    chi1::randType

    tab::tabType

    gtmp::rateNoiseType
    ktmp::rateType

    Y100::uType
    Y200::uType
    Y300::uType
    Y400::uType
    Y1jajb::Array{uType}
    Y2jajb::Array{uType}
    Y3jajb::Array{uType}
    Y4jajb::Array{uType}

    tmpu::uType
    #tmpg::rateNoiseType

end

function alg_cache(
        alg::NON, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = NONConstantCache(real(uBottomEltypeNoUnits))

    gtmp = zero(noise_rate_prototype)
    ktmp = zero(rate_prototype)

    Y100 = zero(u)
    Y200 = zero(u)
    Y300 = zero(u)
    Y400 = zero(u)
    Y1jajb = [zero(u) for ja in 1:m, jb in 1:m]
    Y2jajb = [zero(u) for ja in 1:m, jb in 1:m]
    Y3jajb = [zero(u) for ja in 1:m, jb in 1:m]
    Y4jajb = [zero(u) for ja in 1:m, jb in 1:m]

    tmpu = zero(u)

    return NONCache(
        u, uprev, _dW, _dZ, chi1, tab, gtmp, ktmp, Y100, Y200,
        Y300, Y400, Y1jajb, Y2jajb, Y3jajb, Y4jajb, tmpu
    )
end

# COM
struct COMConstantCache{T} <: StochasticDiffEqConstantCache
    c01::T
    c02::T
    c03::T
    c04::T

    cj1::T
    cj2::T
    cj3::T
    cj4::T

    a0021::T
    a0032::T
    a0043::T

    aj021::T
    aj041::T

    a0j21::T
    a0j31::T
    a0j32::T
    a0j41::T

    ajj21::T
    ajj31::T
    ajj32::T
    ajj41::T
    ajj42::T
    ajj43::T

    ajl31::T
    ajl32::T
    ajl41::T
    ajl42::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function COMConstantCache(T::Type)
    c01 = convert(T, 1 // 6)
    c02 = convert(T, 1 // 3)
    c03 = convert(T, 1 // 3)
    c04 = convert(T, 1 // 6)

    cj1 = convert(T, 1 // 8)
    cj2 = convert(T, 3 // 8)
    cj3 = convert(T, 3 // 8)
    cj4 = convert(T, 1 // 8)

    a0021 = convert(T, 1 // 2)
    a0032 = convert(T, 1 // 2)
    a0043 = convert(T, 1)

    aj021 = convert(T, 2)
    aj041 = convert(T, -2)

    a0j21 = convert(T, 1)
    a0j31 = convert(T, -9 // 8)
    a0j32 = convert(T, 9 // 8)
    a0j41 = convert(T, 1)

    ajj21 = convert(T, 2 // 3)
    ajj31 = convert(T, 1 // 12)
    ajj32 = convert(T, 1 // 4)
    ajj41 = convert(T, -5 // 4)
    ajj42 = convert(T, 1 // 4)
    ajj43 = convert(T, 2)

    ajl31 = convert(T, 1 // 4)
    ajl32 = convert(T, 3 // 4)
    ajl41 = convert(T, 1 // 4)
    ajl42 = convert(T, 3 // 4)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return COMConstantCache(
        c01, c02, c03, c04, cj1, cj2, cj3, cj4, a0021, a0032, a0043, aj021,
        aj041, a0j21, a0j31, a0j32, a0j41, ajj21, ajj31, ajj32, ajj41,
        ajj42, ajj43, ajl31, ajl32, ajl41, ajl42, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::COM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return COMConstantCache(real(uBottomEltypeNoUnits))
end

@cache struct COMCache{uType, randType, tabType, rateNoiseType, rateType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType

    tab::tabType

    gtmp::rateNoiseType
    ktmp::rateType

    Y10::uType
    Y20::uType
    Y30::uType
    Y40::uType
    Y1j::rateNoiseType
    Y2j::rateNoiseType
    Y3j::rateNoiseType
    Y4j::rateNoiseType

    tmpu::uType
    tmpu2::uType
end

function alg_cache(
        alg::COM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
    else
        _dW = zero(ΔW)
    end
    m = length(ΔW)
    tab = COMConstantCache(real(uBottomEltypeNoUnits))

    gtmp = zero(noise_rate_prototype)
    ktmp = zero(rate_prototype)

    Y10 = zero(u)
    Y20 = zero(u)
    Y30 = zero(u)
    Y40 = zero(u)
    Y1j = zero(noise_rate_prototype)
    Y2j = zero(noise_rate_prototype)
    Y3j = zero(noise_rate_prototype)
    Y4j = zero(noise_rate_prototype)

    tmpu = zero(u)
    tmpu2 = zero(u)

    return COMCache(
        u, uprev, _dW, tab, gtmp, ktmp, Y10, Y20, Y30, Y40, Y1j, Y2j, Y3j, Y4j, tmpu, tmpu2
    )
end

# NON2
struct NON2ConstantCache{T} <: StochasticDiffEqConstantCache
    c01::T
    c02::T
    c03::T
    c04::T

    cj1::T
    cj2::T
    cj3::T
    cj4::T

    a0021::T
    a0032::T
    a0043::T

    aj021::T
    aj041::T

    a0j21::T
    a0j31::T
    a0j32::T
    a0j41::T

    ajj21::T
    ajj31::T
    ajj32::T
    ajj41::T
    ajj42::T
    ajj43::T

    ajl31::T
    ajl32::T
    ajl41::T
    ajl42::T

    # for non-commuting terms
    #ckj3::T
    #ckj4::T
    #akjjl32::T
    #akjjl42::T
    # are all expressed in terms of γ
    γ::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function NON2ConstantCache(T::Type)
    c01 = convert(T, 1 // 6)
    c02 = convert(T, 1 // 3)
    c03 = convert(T, 1 // 3)
    c04 = convert(T, 1 // 6)

    cj1 = convert(T, 1 // 8)
    cj2 = convert(T, 3 // 8)
    cj3 = convert(T, 3 // 8)
    cj4 = convert(T, 1 // 8)

    a0021 = convert(T, 1 // 2)
    a0032 = convert(T, 1 // 2)
    a0043 = convert(T, 1)

    aj021 = convert(T, 2)
    aj041 = convert(T, -2)

    a0j21 = convert(T, 1)
    a0j31 = convert(T, -9 // 8)
    a0j32 = convert(T, 9 // 8)
    a0j41 = convert(T, 1)

    ajj21 = convert(T, 2 // 3)
    ajj31 = convert(T, 1 // 12)
    ajj32 = convert(T, 1 // 4)
    ajj41 = convert(T, -5 // 4)
    ajj42 = convert(T, 1 // 4)
    ajj43 = convert(T, 2)

    ajl31 = convert(T, 1 // 4)
    ajl32 = convert(T, 3 // 4)
    ajl41 = convert(T, 1 // 4)
    ajl42 = convert(T, 3 // 4)

    γ = convert(T, 1)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return NON2ConstantCache(
        c01, c02, c03, c04, cj1, cj2, cj3, cj4, a0021, a0032, a0043, aj021,
        aj041, a0j21, a0j31, a0j32, a0j41, ajj21, ajj31, ajj32, ajj41,
        ajj42, ajj43, ajl31, ajl32, ajl41, ajl42, γ, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::NON2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NON2ConstantCache(real(uBottomEltypeNoUnits))
end

@cache struct NON2Cache{uType, randType, MType1, tabType, rateNoiseType, rateType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    _dW::randType
    _dZ::randType
    chi1::randType
    Ihat2::MType1

    tab::tabType

    gtmp::rateNoiseType
    gtmp1::rateNoiseType
    ktmp::rateType

    Y10::uType
    Y20::uType
    Y30::uType
    Y40::uType
    Y1j::rateNoiseType
    Y2j::rateNoiseType
    Y3j::rateNoiseType
    Y4j::rateNoiseType

    tmpu::uType
    tmpu2::uType
end

function alg_cache(
        alg::NON2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = copy(ΔW)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zero(ΔW)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    Ihat2 = zeros(eltype(ΔW), m, m)
    tab = NON2ConstantCache(real(uBottomEltypeNoUnits))

    gtmp = zero(noise_rate_prototype)
    gtmp1 = zero(noise_rate_prototype)
    ktmp = zero(rate_prototype)

    Y10 = zero(u)
    Y20 = zero(u)
    Y30 = zero(u)
    Y40 = zero(u)
    Y1j = zero(noise_rate_prototype)
    Y2j = zero(noise_rate_prototype)
    Y3j = zero(noise_rate_prototype)
    Y4j = zero(noise_rate_prototype)

    tmpu = zero(u)
    tmpu2 = zero(u)

    return NON2Cache(
        u, uprev, _dW, _dZ, chi1, Ihat2, tab, gtmp, gtmp1, ktmp,
        Y10, Y20, Y30, Y40, Y1j, Y2j, Y3j, Y4j, tmpu, tmpu2
    )
end

# SIE / SME methods

struct SIESMEConstantCache{T, T2} <: StochasticDiffEqConstantCache
    α1::T
    α2::T

    γ1::T

    λ1::T
    λ2::T
    λ3::T

    µ1::T
    µ2::T
    µ3::T

    µ0::T2
    µbar0::T2

    λ0::T
    λbar0::T

    ν1::T
    ν2::T

    β2::T
    β3::T

    δ2::T
    δ3::T
end

function SIEAConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    α1 = convert(T, 1 // 2)
    α2 = convert(T, 1 // 2)

    γ1 = convert(T, 1 // 2)

    λ1 = convert(T, 1 // 4)
    λ2 = convert(T, -1 // 4)
    λ3 = convert(T, 1 // 4)

    µ1 = convert(T, 1 // 4)
    µ2 = convert(T, 1 // 4)
    µ3 = convert(T, -1 // 4)

    µ0 = convert(T2, 1 // 1)
    µbar0 = convert(T2, 1 // 1)

    λ0 = convert(T, 1 // 1)
    λbar0 = convert(T, 1 // 1)

    ν1 = convert(T, 1 // 1)
    ν2 = convert(T, 0)

    β2 = convert(T, 1 // 1)
    β3 = convert(T, 0)

    δ2 = convert(T, -1 // 1)
    δ3 = convert(T, 0)

    return SIESMEConstantCache(
        α1, α2, γ1, λ1, λ2, λ3, µ1, µ2, µ3, µ0, µbar0, λ0, λbar0, ν1, ν2, β2, β3, δ2, δ3
    )
end

function alg_cache(
        alg::SIEA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SIEAConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function SMEAConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    α1 = convert(T, 0)
    α2 = convert(T, 1 // 1)

    γ1 = convert(T, 1 // 2)

    λ1 = convert(T, 1 // 4)
    λ2 = convert(T, -1 // 4)
    λ3 = convert(T, 1 // 4)

    µ1 = convert(T, 1 // 4)
    µ2 = convert(T, 1 // 4)
    µ3 = convert(T, -1 // 4)

    µ0 = convert(T2, 1 // 2)
    µbar0 = convert(T2, 1 // 1)

    λ0 = convert(T, 1 // 2)
    λbar0 = convert(T, 1 // 1)

    ν1 = convert(T, (2 - sqrt(6)) / 4)
    ν2 = convert(T, sqrt(6) / 12)

    β2 = convert(T, 1 // 1)
    β3 = convert(T, 0)

    δ2 = convert(T, -1 // 1)
    δ3 = convert(T, 0)

    return SIESMEConstantCache(
        α1, α2, γ1, λ1, λ2, λ3, µ1, µ2, µ3, µ0, µbar0, λ0, λbar0, ν1, ν2, β2, β3, δ2, δ3
    )
end

function alg_cache(
        alg::SMEA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SMEAConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function SIEBConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    α1 = convert(T, 1 // 2)
    α2 = convert(T, 1 // 2)

    γ1 = convert(T, -1 // 5)

    λ1 = convert(T, 3 // 5)
    λ2 = convert(T, 3 // 2)
    λ3 = convert(T, -1 // 2)

    µ1 = convert(T, 3 // 5)
    µ2 = convert(T, -3 // 2)
    µ3 = convert(T, 1 // 2)

    µ0 = convert(T2, 1 // 1)
    µbar0 = convert(T2, 5 // 12)

    λ0 = convert(T, 1 // 1)
    λbar0 = convert(T, 5 // 12)

    ν1 = convert(T, 1 // 1)
    ν2 = convert(T, 0)

    β2 = convert(T, 0)
    β3 = convert(T, -1 // 6)

    δ2 = convert(T, 0)
    δ3 = convert(T, 1 // 6)

    return SIESMEConstantCache(
        α1, α2, γ1, λ1, λ2, λ3, µ1, µ2, µ3, µ0, µbar0, λ0, λbar0, ν1, ν2, β2, β3, δ2, δ3
    )
end

function alg_cache(
        alg::SIEB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SIEBConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function SMEBConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    α1 = convert(T, 0)
    α2 = convert(T, 1 // 1)

    γ1 = convert(T, -1 // 5)

    λ1 = convert(T, 3 // 5)
    λ2 = convert(T, 3 // 2)
    λ3 = convert(T, -1 // 2)

    µ1 = convert(T, 3 // 5)
    µ2 = convert(T, -3 // 2)
    µ3 = convert(T, 1 // 2)

    µ0 = convert(T2, 1 // 2)
    µbar0 = convert(T2, 5 // 12)

    λ0 = convert(T, 1 // 2)
    λbar0 = convert(T, 5 // 12)

    ν1 = convert(T, (2 - sqrt(6)) / 4)
    ν2 = convert(T, sqrt(6) / 12)

    β2 = convert(T, 0)
    β3 = convert(T, -1 // 6)

    δ2 = convert(T, 0)
    δ3 = convert(T, 1 // 6)

    return SIESMEConstantCache(
        α1, α2, γ1, λ1, λ2, λ3, µ1, µ2, µ3, µ0, µbar0, λ0, λbar0, ν1, ν2, β2, β3, δ2, δ3
    )
end

function alg_cache(
        alg::SMEB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SMEBConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SIESMECache{uType, randType, tabType, rateNoiseType, rateType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType

    W2::randType
    W3::randType

    tab::tabType

    k0::rateType
    k1::rateType

    g0::rateNoiseType
    g1::rateNoiseType
    g2::rateNoiseType

    #tmp1::possibleRateType
    tmpu::uType
end

function alg_cache(
        alg::SIEA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        W2 = copy(ΔW)
        W3 = copy(ΔW)
    else
        W2 = zero(ΔW)
        W3 = zero(ΔW)
    end

    k0 = zero(u)
    k1 = zero(u)

    g0 = zero(noise_rate_prototype)
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)

    tmpu = zero(u)

    tab = SIEAConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))

    return SIESMECache(u, uprev, W2, W3, tab, k0, k1, g0, g1, g2, tmpu)
end

function alg_cache(
        alg::SMEA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        W2 = copy(ΔW)
        W3 = copy(ΔW)
    else
        W2 = zero(ΔW)
        W3 = zero(ΔW)
    end

    k0 = zero(u)
    k1 = zero(u)

    g0 = zero(noise_rate_prototype)
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)

    tmpu = zero(u)

    tab = SMEAConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))

    return SIESMECache(u, uprev, W2, W3, tab, k0, k1, g0, g1, g2, tmpu)
end

function alg_cache(
        alg::SIEB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        W2 = copy(ΔW)
        W3 = copy(ΔW)
    else
        W2 = zero(ΔW)
        W3 = zero(ΔW)
    end

    k0 = zero(u)
    k1 = zero(u)

    g0 = zero(noise_rate_prototype)
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)

    tmpu = zero(u)

    tab = SIEBConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))

    return SIESMECache(u, uprev, W2, W3, tab, k0, k1, g0, g1, g2, tmpu)
end

function alg_cache(
        alg::SMEB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        W2 = copy(ΔW)
        W3 = copy(ΔW)
    else
        W2 = zero(ΔW)
        W3 = zero(ΔW)
    end

    k0 = zero(u)
    k1 = zero(u)

    g0 = zero(noise_rate_prototype)
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)

    tmpu = zero(u)

    tab = SMEBConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))

    return SIESMECache(u, uprev, W2, W3, tab, k0, k1, g0, g1, g2, tmpu)
end

# Tang & Xiao: DOI 10.1007/s10543-016-0618-9 W2Ito1 and W2Ito2 methods

struct W2Ito1ConstantCache{T} <: StochasticDiffEqConstantCache
    # hard-coded version
    a021::T
    a031::T
    a032::T

    a121::T
    a131::T
    #a132::T

    #a221::T
    #a231::T
    #a232::T

    b021::T
    b031::T
    #b032::T

    b121::T
    b131::T
    #b132::T

    b221::T
    #b222::T
    #b223::T
    #b231::T
    #b232::T
    #b233::T

    α1::T
    α2::T
    α3::T

    beta01::T
    beta02::T
    beta03::T

    beta11::T
    #beta12::T
    beta13::T

    #quantile(Normal(),1/6)
    NORMAL_ONESIX_QUANTILE::T
end

function W2Ito1ConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
    a021 = convert(T, 1 // 2)
    a031 = convert(T, -1)
    a032 = convert(T, 2)

    a121 = convert(T, 1 // 4)
    a131 = convert(T, 1 // 4)

    b021 = convert(T, (6 - sqrt(6)) / 10)
    b031 = convert(T, (3 + 2 * sqrt(6)) / 5)

    b121 = convert(T, 1 // 2)
    b131 = convert(T, -1 // 2)

    b221 = convert(T, 1)

    α1 = convert(T, 1 // 6)
    α2 = convert(T, 2 // 3)
    α3 = convert(T, 1 // 6)

    beta01 = convert(T, -1)
    beta02 = convert(T, 1)
    beta03 = convert(T, 1)

    beta11 = convert(T, 2)
    beta13 = convert(T, -2)

    NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

    return W2Ito1ConstantCache(
        a021, a031, a032, a121, a131, b021, b031, b121, b131, b221, α1, α2,
        α3, beta01, beta02, beta03, beta11, beta13, NORMAL_ONESIX_QUANTILE
    )
end

function alg_cache(
        alg::W2Ito1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return W2Ito1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct W2Ito1Cache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhat::uType

    _dW::randType
    _dZ::randType
    chi1::randType

    tab::tabType

    g1::rateNoiseType
    g2::rateNoiseType
    g3::rateNoiseType

    k1::rateType
    k2::rateType
    k3::rateType

    H02::uType
    H03::uType
    H12::Vector{uType}
    H13::Vector{uType}

    tmp1::possibleRateType
    tmpg::rateNoiseType

    tmp::uType
    resids::uType
end

function alg_cache(
        alg::W2Ito1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        _dW = copy(ΔW)
        _dZ = zeros(eltype(ΔW), 2)
        chi1 = copy(ΔW)
    else
        _dW = zero(ΔW)
        _dZ = zeros(eltype(ΔW), 2)
        chi1 = zero(ΔW)
    end
    m = length(ΔW)
    tab = W2Ito1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)

    H02 = zero(u)
    H03 = zero(u)
    H12 = Vector{typeof(u)}()
    H13 = Vector{typeof(u)}()

    for k in 1:m
        push!(H12, zero(u))
        push!(H13, zero(u))
    end

    tmp1 = zero(rate_prototype)
    tmpg = zero(noise_rate_prototype)

    uhat = copy(uprev)
    tmp = zero(u)
    resids = zero(u)

    return W2Ito1Cache(
        u, uprev, uhat, _dW, _dZ, chi1, tab, g1, g2, g3, k1,
        k2, k3, H02, H03, H12, H13, tmp1, tmpg, tmp, resids
    )
end
