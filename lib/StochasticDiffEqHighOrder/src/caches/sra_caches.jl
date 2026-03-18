struct SRA1ConstantCache <: StochasticDiffEqConstantCache end
function alg_cache(
        alg::SRA1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRA1ConstantCache()
end

@cache struct SRA1Cache{randType, rateType, uType, rateNoiseType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    chi2::randType
    tmp1::uType
    E₁::rateType
    E₂::rateType
    gt::rateNoiseType
    k₁::rateType
    k₂::rateType
    gpdt::rateNoiseType
    tmp::uType
end

function alg_cache(
        alg::SRA1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi2 = copy(ΔW)
    else
        chi2 = zero(ΔW)
    end
    tmp1 = zero(u)
    E₁ = zero(rate_prototype)
    gt = zero(noise_rate_prototype)
    gpdt = zero(noise_rate_prototype)
    E₂ = zero(rate_prototype)
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    tmp = zero(u)
    return SRA1Cache(u, uprev, chi2, tmp1, E₁, E₂, gt, k₁, k₂, gpdt, tmp)
end

struct SRA2ConstantCache{T, T2} <: StochasticDiffEqConstantCache
    a21::T
    b21::T
    c02::T2
    c11::T2
    c12::T2
    α1::T
    α2::T
    beta12::T
    beta21::T
    beta22::T
end

function SRA2ConstantCache(T::Type, T2::Type)
    a21 = convert(T, 3 // 4)
    b21 = convert(T, 3 // 2)
    c02 = convert(T2, 3 // 4)
    c11 = convert(T2, 1 // 3)
    c12 = convert(T2, 1)
    α1 = convert(T, 1 // 3)
    α2 = convert(T, 2 // 3)
    beta12 = convert(T, 1)
    beta21 = convert(T, 3 // 2)
    beta22 = convert(T, -3 // 2)
    return SRA2ConstantCache(a21, b21, c02, c11, c12, α1, α2, beta12, beta21, beta22)
end

function alg_cache(
        alg::SRA2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRA2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SRA2Cache{uType, randType, tabType, rateNoiseType, T} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    chi2::randType
    tab::tabType
    g1::rateNoiseType
    g2::rateNoiseType
    k1::T
    k2::T
    E₁::T
    E₂::T
    tmp::T
end

function alg_cache(
        alg::SRA2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi2 = copy(ΔW)
    else
        chi2 = zero(ΔW)
    end
    tab = SRA2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = k1
    return SRA2Cache(u, uprev, chi2, tab, g1, g2, k1, k2, E₁, E₂, tmp)
end

struct ThreeStageSRAConstantCache{T, T2} <: StochasticDiffEqConstantCache
    a21::T
    a31::T
    a32::T
    b21::T
    b31::T
    b32::T
    c02::T2
    c03::T2
    c11::T2
    c12::T2
    c13::T2
    α1::T
    α2::T
    α3::T
    beta11::T
    beta12::T
    beta13::T
    beta21::T
    beta22::T
    beta23::T
end

function SRA3ConstantCache(T::Type, T2::Type)
    a21 = convert(T, 1)
    a31 = convert(T, 1 // 4)
    a32 = convert(T, 1 // 4)
    b21 = convert(T, 0)
    b31 = convert(T, 1)
    b32 = convert(T, 1 // 2)
    c02 = convert(T2, 1)
    c03 = convert(T2, 1 // 2)
    c11 = convert(T2, 1)
    c12 = convert(T2, 0)
    c13 = convert(T2, 0)
    α1 = convert(T, 1 // 6)
    α2 = convert(T, 1 // 6)
    α3 = convert(T, 2 // 3)
    beta11 = convert(T, 1)
    beta12 = convert(T, 0)
    beta13 = convert(T, 0)
    beta21 = convert(T, -1)
    beta22 = convert(T, 1)
    beta23 = convert(T, 0)
    return ThreeStageSRAConstantCache(
        a21, a31, a32, b21, b31, b32, c02, c03, c11, c12, c13,
        α1, α2, α3, beta11, beta12, beta13, beta21, beta22, beta23
    )
end

function SOSRAConstantCache(T::Type, T2::Type)
    α1 = convert(T, 0.2889874966892885)
    α2 = convert(T, 0.6859880440839937)
    α3 = convert(T, 0.025024459226717772)
    c02 = convert(T2, 0.6923962376159507)
    c03 = convert(T2, 1)
    c11 = convert(T2, 0)
    c12 = convert(T2, 0.041248171110700504)
    c13 = convert(T2, 1)
    beta11 = convert(T, -16.792534242221663)
    beta12 = convert(T, 17.514995785380226)
    beta13 = convert(T, 0.27753845684143835)
    beta21 = convert(T, 0.4237535769069274)
    beta22 = convert(T, 0.6010381474428539)
    beta23 = convert(T, -1.0247917243497813)
    a21 = convert(T, 0.6923962376159507)
    a31 = convert(T, -3.1609142252828395)
    a32 = convert(T, 4.1609142252828395)
    b21 = convert(T, 1.3371632704399763)
    b31 = convert(T, 1.442371048468624)
    b32 = convert(T, 1.8632741501139225)
    return ThreeStageSRAConstantCache(
        a21, a31, a32, b21, b31, b32, c02, c03, c11, c12, c13,
        α1, α2, α3, beta11, beta12, beta13, beta21, beta22, beta23
    )
end

function SOSRA2ConstantCache(T::Type, T2::Type)
    α1 = convert(T, 0.4999999999999998)
    α2 = convert(T, -0.9683897375354181)
    α3 = convert(T, 1.4683897375354185)
    c02 = convert(T2, 1)
    c03 = convert(T2, 1)
    c11 = convert(T2, 0)
    c12 = convert(T2, 1)
    c13 = convert(T2, 1)
    beta11 = convert(T, 0)
    beta12 = convert(T, 0.92438032145683)
    beta13 = convert(T, 0.07561967854316998)
    beta21 = convert(T, 1)
    beta22 = convert(T, -0.8169981105823436)
    beta23 = convert(T, -0.18300188941765633)
    a21 = convert(T, 1)
    a31 = convert(T, 0.9511849235504364)
    a32 = convert(T, 0.04881507644956362)
    b21 = convert(T, 0.7686101171003622)
    b31 = convert(T, 0.43886792994934987)
    b32 = convert(T, 0.7490415909204886)
    return ThreeStageSRAConstantCache(
        a21, a31, a32, b21, b31, b32, c02, c03, c11, c12, c13,
        α1, α2, α3, beta11, beta12, beta13, beta21, beta22, beta23
    )
end

function alg_cache(
        alg::SRA3, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
        f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRA3ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function alg_cache(
        alg::SOSRA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
        f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SOSRAConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function alg_cache(
        alg::SOSRA2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
        f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SOSRA2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct ThreeStageSRACache{uType, randType, tabType, rateNoiseType, rateType, GT} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    chi2::randType
    tab::tabType
    g1::rateNoiseType
    g2::rateNoiseType
    g3::rateNoiseType
    k1::rateType
    k2::rateType
    k3::rateType
    E₁::rateType
    E₂::rateType
    tmp::rateType
    gtmp::GT
end

function alg_cache(
        alg::SRA3, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi2 = copy(ΔW)
    else
        chi2 = zero(ΔW)
    end
    tab = SRA3ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = k1
    if typeof(noise_rate_prototype) == typeof(rate_prototype)
        gtmp = nothing
    else
        gtmp = zero(noise_rate_prototype)
    end
    return ThreeStageSRACache(u, uprev, chi2, tab, g1, g2, g3, k1, k2, k3, E₁, E₂, tmp, gtmp)
end

function alg_cache(
        alg::SOSRA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi2 = copy(ΔW)
    else
        chi2 = zero(ΔW)
    end
    tab = SOSRAConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = k1
    if typeof(noise_rate_prototype) == typeof(rate_prototype)
        gtmp = nothing
    else
        gtmp = zero(noise_rate_prototype)
    end
    return ThreeStageSRACache(u, uprev, chi2, tab, g1, g2, g3, k1, k2, k3, E₁, E₂, tmp, gtmp)
end

function alg_cache(
        alg::SOSRA2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi2 = copy(ΔW)
    else
        chi2 = zero(ΔW)
    end
    tab = SOSRA2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = k1
    if typeof(noise_rate_prototype) == typeof(rate_prototype)
        gtmp = nothing
    else
        gtmp = zero(noise_rate_prototype)
    end
    return ThreeStageSRACache(u, uprev, chi2, tab, g1, g2, g3, k1, k2, k3, E₁, E₂, tmp, gtmp)
end

struct SRAConstantCache{VType1, VType2, MType, uType} <: StochasticDiffEqConstantCache
    c₀::VType1
    c₁::VType1
    A₀::MType
    B₀::MType
    α::VType2
    β₁::VType2
    β₂::VType2
    stages::Int
    H0::Vector{uType}
end

function SRAConstantCache(tableau, rate_prototype)
    (; c₀, c₁, A₀, B₀, α, β₁, β₂) = tableau
    stages = length(α)
    H0 = Vector{typeof(rate_prototype)}(undef, stages)
    return SRAConstantCache(c₀, c₁, A₀', B₀', α, β₁, β₂, stages, H0)
end

function alg_cache(
        alg::SRA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRAConstantCache(alg.tableau, rate_prototype)
end

@cache struct SRACache{uType, rateType, tabType, randType, rateNoiseType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    H0::Vector{uType}
    A0temp::rateType
    B0temp::rateType
    ftmp::rateType
    gtmp::rateNoiseType
    chi2::randType
    atemp::rateType
    btemp::rateType
    E₁::rateType
    E₁temp::rateType
    E₂::rateType
    tmp::uType
    tab::tabType
end

u_cache(c::SRACache) = ()
function du_cache(c::SRACache)
    return (
        c.A0temp, c.B0temp, c.ftmp, c.gtmp, c.chi2, c.chi2, c.atemp,
        c.btemp, c.E₁, c.E₁temp, c.E₂,
    )
end
user_cache(c::SRACache) = (c.u, c.uprev, c.tmp, c.H0...)

function alg_cache(
        alg::SRA, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    H0 = Vector{typeof(u)}()
    tab = SRAConstantCache(alg.tableau, rate_prototype)
    for i in 1:tab.stages
        push!(H0, zero(u))
    end
    A0temp = zero(rate_prototype)
    B0temp = zero(rate_prototype)
    ftmp = zero(rate_prototype)
    gtmp = zero(noise_rate_prototype)
    if ΔW isa Union{SArray, Number}
        chi2 = copy(ΔW)
    else
        chi2 = zero(ΔW)
    end
    atemp = zero(rate_prototype)
    btemp = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    E₁temp = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    tmp = zero(u)
    return SRACache(
        u, uprev, H0, A0temp, B0temp, ftmp, gtmp,
        chi2, atemp, btemp, E₁, E₁temp, E₂, tmp, tab
    )
end
