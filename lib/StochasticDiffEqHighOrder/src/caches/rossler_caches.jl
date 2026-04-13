struct SRIConstantCache{VType1, VType2, MType, uType} <: StochasticDiffEqConstantCache
    c₀::VType1
    c₁::VType1
    A₀::MType
    A₁::MType
    B₀::MType
    B₁::MType
    α::VType2
    β₁::VType2
    β₂::VType2
    β₃::VType2
    β₄::VType2
    stages::Int
    H0::Vector{uType}
    H1::Vector{uType}
    error_terms::Int
end

function SRIConstantCache(tableau, rate_prototype, error_terms)
    (; c₀, c₁, A₀, A₁, B₀, B₁, α, β₁, β₂, β₃, β₄) = tableau
    stages = length(α)
    H0 = Array{typeof(rate_prototype)}(undef, stages)
    H1 = Array{typeof(rate_prototype)}(undef, stages)
    return SRIConstantCache(
        c₀, c₁, A₀', A₁', B₀', B₁', α, β₁, β₂, β₃, β₄, stages, H0, H1, error_terms
    )
end

function alg_cache(
        alg::SRI, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRIConstantCache(alg.tableau, rate_prototype, alg.error_terms)
end

@cache struct SRICache{randType, uType, rateType, tabType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    H0::Vector{uType}
    H1::Vector{uType}
    A0temp::rateType
    A1temp::rateType
    B0temp::rateType
    B1temp::rateType
    A0temp2::rateType
    A1temp2::rateType
    B0temp2::rateType
    B1temp2::rateType
    atemp::rateType
    btemp::rateType
    E₁::rateType
    E₂::rateType
    E₁temp::rateType
    ftemp::rateType
    gtemp::rateType
    chi1::randType
    chi2::randType
    chi3::randType
    tmp::uType
    tab::tabType
end

function alg_cache(
        alg::SRI, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    H0 = Vector{typeof(u)}()
    H1 = Vector{typeof(u)}()
    tab = SRIConstantCache(alg.tableau, rate_prototype, alg.error_terms)
    for i in 1:tab.stages
        push!(H0, zero(u))
        push!(H1, zero(u))
    end
    #TODO Reduce memory
    A0temp = zero(rate_prototype)
    A1temp = zero(rate_prototype)
    B0temp = zero(rate_prototype)
    B1temp = zero(rate_prototype)
    A0temp2 = zero(rate_prototype)
    A1temp2 = zero(rate_prototype)
    B0temp2 = zero(rate_prototype)
    B1temp2 = zero(rate_prototype)
    atemp = zero(rate_prototype)
    btemp = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    E₁temp = zero(rate_prototype)
    ftemp = zero(rate_prototype)
    gtemp = zero(rate_prototype)

    if ΔW isa Union{SArray, Number}
        chi1 = copy(ΔW)
        chi2 = copy(ΔW)
        chi3 = copy(ΔW)
    else
        chi1 = zero(ΔW)
        chi2 = zero(ΔW)
        chi3 = zero(ΔW)
    end
    tmp = zero(u)
    return SRICache(
        u, uprev, H0, H1, A0temp, A1temp, B0temp,
        B1temp, A0temp2, A1temp2, B0temp2, B1temp2,
        atemp, btemp, E₁, E₂, E₁temp, ftemp, gtemp, chi1, chi2, chi3, tmp, tab
    )
end

struct SRIW1ConstantCache <: StochasticDiffEqConstantCache end
function alg_cache(
        alg::SRIW1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRIW1ConstantCache()
end

@cache struct SRIW1Cache{randType, uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    chi1::randType
    chi2::randType
    chi3::randType
    fH01o4::rateType
    g₁o2::rateType
    H0::uType
    H11::uType
    H12::uType
    H13::uType
    g₂o3::rateType
    Fg₂o3::rateType
    g₃o3::rateType
    Tg₃o3::rateType
    mg₁::rateType
    E₁::rateType
    E₂::rateType
    fH01::rateType
    fH02::rateType
    g₁::rateType
    g₂::rateType
    g₃::rateType
    g₄::rateType
    tmp::uType
end

function alg_cache(
        alg::SRIW1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi1 = copy(ΔW)
        chi2 = copy(ΔW)
        chi3 = copy(ΔW)
    else
        chi1 = zero(ΔW)
        chi2 = zero(ΔW)
        chi3 = zero(ΔW)
    end
    fH01o4 = zero(rate_prototype)
    g₁o2 = zero(rate_prototype)
    H0 = zero(u)
    H11 = zero(u)
    H12 = zero(u)
    H13 = zero(u)
    g₂o3 = zero(rate_prototype)
    Fg₂o3 = zero(rate_prototype)
    g₃o3 = zero(rate_prototype)
    Tg₃o3 = zero(rate_prototype)
    mg₁ = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    fH01 = zero(rate_prototype)
    fH02 = zero(rate_prototype)
    g₁ = zero(rate_prototype)
    g₂ = zero(rate_prototype)
    g₃ = zero(rate_prototype)
    g₄ = zero(rate_prototype)
    tmp = zero(u)
    return SRIW1Cache(
        u, uprev, chi1, chi2, chi3, fH01o4, g₁o2, H0, H11, H12, H13, g₂o3,
        Fg₂o3, g₃o3, Tg₃o3, mg₁, E₁, E₂, fH01, fH02, g₁, g₂, g₃, g₄, tmp
    )
end

struct FourStageSRIConstantCache{T, T2} <: StochasticDiffEqConstantCache
    a021::T
    a031::T
    a032::T
    a041::T
    a042::T
    a043::T
    a121::T
    a131::T
    a132::T
    a141::T
    a142::T
    a143::T
    b021::T
    b031::T
    b032::T
    b041::T
    b042::T
    b043::T
    b121::T
    b131::T
    b132::T
    b141::T
    b142::T
    b143::T
    α1::T
    α2::T
    α3::T
    α4::T
    c02::T2
    c03::T2
    c04::T2
    c11::T2
    c12::T2
    c13::T2
    c14::T2
    beta11::T
    beta12::T
    beta13::T
    beta14::T
    beta21::T
    beta22::T
    beta23::T
    beta24::T
    beta31::T
    beta32::T
    beta33::T
    beta34::T
    beta41::T
    beta42::T
    beta43::T
    beta44::T
end

function SRIW2ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 1)
    a031 = convert(T, 1 // 4)
    a032 = convert(T, 1 // 4)
    a041 = convert(T, 0)
    a042 = convert(T, 0)
    a043 = convert(T, 0)
    a121 = convert(T, 1 // 4)
    a131 = convert(T, 1)
    a132 = convert(T, 0)
    a141 = convert(T, 0)
    a142 = convert(T, 0)
    a143 = convert(T, 1 // 4)
    b021 = convert(T, 0)
    b031 = convert(T, 1)
    b032 = convert(T, 1 // 2)
    b041 = convert(T, 0)
    b042 = convert(T, 0)
    b043 = convert(T, 0)
    b121 = convert(T, -1 // 2)
    b131 = convert(T, 1)
    b132 = convert(T, 0)
    b141 = convert(T, 2)
    b142 = convert(T, -1)
    b143 = convert(T, 1 // 2)
    α1 = convert(T, 1 // 6)
    α2 = convert(T, 1 // 6)
    α3 = convert(T, 2 // 3)
    α4 = convert(T, 0)
    c02 = convert(T2, 1)
    c03 = convert(T2, 1 // 2)
    c04 = convert(T2, 0)
    c11 = convert(T2, 0)
    c12 = convert(T2, 1 // 4)
    c13 = convert(T2, 1)
    c14 = convert(T2, 1 // 4)
    beta11 = convert(T, -1)
    beta12 = convert(T, 4 // 3)
    beta13 = convert(T, 2 // 3)
    beta14 = convert(T, 0)
    beta21 = convert(T, 1)
    beta22 = convert(T, -4 // 3)
    beta23 = convert(T, 1 // 3)
    beta24 = convert(T, 0)
    beta31 = convert(T, 2)
    beta32 = convert(T, -4 // 3)
    beta33 = convert(T, -2 // 3)
    beta34 = convert(T, 0)
    beta41 = convert(T, -2)
    beta42 = convert(T, 5 // 3)
    beta43 = convert(T, -2 // 3)
    beta44 = convert(T, 1)
    return FourStageSRIConstantCache(
        a021, a031, a032, a041, a042, a043, a121, a131, a132, a141,
        a142, a143, b021, b031, b032, b041, b042, b043,
        b121, b131, b132, b141, b142, b143, α1, α2,
        α3, α4, c02, c03, c04, c11, c12, c13, c14, beta11,
        beta12, beta13, beta14, beta21, beta22,
        beta23, beta24, beta31, beta32, beta33, beta34, beta41, beta42, beta43, beta44
    )
end

function alg_cache(
        alg::SRIW2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SRIW2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function SOSRIConstantCache(T::Type, T2::Type)
    a021 = convert(T, -0.04199224421316468)
    a031 = convert(T, 2.842612915017106)
    a032 = convert(T, -2.0527723684000727)
    a041 = convert(T, 4.338237071435815)
    a042 = convert(T, -2.8895936137439793)
    a043 = convert(T, 2.3017575594644466)
    a121 = convert(T, 0.26204282091330466)
    a131 = convert(T, 0.20903646383505375)
    a132 = convert(T, -0.1502377115150361)
    a141 = convert(T, 0.05836595312746999)
    a142 = convert(T, 0.6149440396332373)
    a143 = convert(T, 0.08535117634046772)
    b021 = convert(T, -0.21641093549612528)
    b031 = convert(T, 1.5336352863679572)
    b032 = convert(T, 0.26066223492647056)
    b041 = convert(T, -1.0536037558179159)
    b042 = convert(T, 1.7015284721089472)
    b043 = convert(T, -0.20725685784180017)
    b121 = convert(T, -0.5119011827621657)
    b131 = convert(T, 2.67767339866713)
    b132 = convert(T, -4.9395031322250995)
    b141 = convert(T, 0.15580956238299215)
    b142 = convert(T, 3.2361551006624674)
    b143 = convert(T, -1.4223118283355949)
    α1 = convert(T, 1.140099274172029)
    α2 = convert(T, -0.6401334255743456)
    α3 = convert(T, 0.4736296532772559)
    α4 = convert(T, 0.026404498125060714)
    c02 = convert(T2, -0.04199224421316468)
    c03 = convert(T2, 0.7898405466170333)
    c04 = convert(T2, 3.7504010171562823)
    c11 = convert(T2, 0.0)
    c12 = convert(T2, 0.26204282091330466)
    c13 = convert(T2, 0.05879875232001766)
    c14 = convert(T2, 0.758661169101175)
    beta11 = convert(T, -1.8453464565104432)
    beta12 = convert(T, 2.688764531100726)
    beta13 = convert(T, -0.2523866501071323)
    beta14 = convert(T, 0.40896857551684956)
    beta21 = convert(T, 0.4969658141589478)
    beta22 = convert(T, -0.5771202869753592)
    beta23 = convert(T, -0.12919702470322217)
    beta24 = convert(T, 0.2093514975196336)
    beta31 = convert(T, 2.8453464565104425)
    beta32 = convert(T, -2.688764531100725)
    beta33 = convert(T, 0.2523866501071322)
    beta34 = convert(T, -0.40896857551684945)
    beta41 = convert(T, 0.11522663875443433)
    beta42 = convert(T, -0.57877086147738)
    beta43 = convert(T, 0.2857851028163886)
    beta44 = convert(T, 0.17775911990655704)
    return FourStageSRIConstantCache(
        a021, a031, a032, a041, a042, a043, a121, a131, a132, a141,
        a142, a143, b021, b031, b032, b041, b042, b043,
        b121, b131, b132, b141, b142, b143, α1, α2,
        α3, α4, c02, c03, c04, c11, c12, c13, c14, beta11,
        beta12, beta13, beta14, beta21, beta22,
        beta23, beta24, beta31, beta32, beta33, beta34, beta41, beta42, beta43, beta44
    )
end

function alg_cache(
        alg::SOSRI, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SOSRIConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function SOSRI2ConstantCache(T::Type, T2::Type)
    a021 = convert(T, 0.13804532298278663)
    a031 = convert(T, 0.5818361298250374)
    a032 = convert(T, 0.4181638701749618)
    a041 = convert(T, 0.4670018408674211)
    a042 = convert(T, 0.8046204792187386)
    a043 = convert(T, -0.27162232008616016)
    a121 = convert(T, 0.45605532163856893)
    a131 = convert(T, 0.7555807846451692)
    a132 = convert(T, 0.24441921535482677)
    a141 = convert(T, 0.6981181143266059)
    a142 = convert(T, 0.3453277086024727)
    a143 = convert(T, -0.04344582292908241)
    b021 = convert(T, 0.08852381537667678)
    b031 = convert(T, 1.0317752458971061)
    b032 = convert(T, 0.4563552922077882)
    b041 = convert(T, 1.73078280444124)
    b042 = convert(T, -0.46089678470929774)
    b043 = convert(T, -0.9637509618944188)
    b121 = convert(T, 0.6753186815412179)
    b131 = convert(T, -0.07452812525785148)
    b132 = convert(T, -0.49783736486149366)
    b141 = convert(T, -0.5591906709928903)
    b142 = convert(T, 0.022696571806569924)
    b143 = convert(T, -0.8984927888368557)
    α1 = convert(T, -0.15036858140642623)
    α2 = convert(T, 0.7545275856696072)
    α3 = convert(T, 0.686995463807979)
    α4 = convert(T, -0.2911544680711602)
    c02 = convert(T2, 0.13804532298278663)
    c03 = convert(T2, 0.9999999999999992)
    c04 = convert(T2, 0.9999999999999994)
    c11 = convert(T2, 0.0)
    c12 = convert(T2, 0.45605532163856893)
    c13 = convert(T2, 0.999999999999996)
    c14 = convert(T2, 0.9999999999999962)
    beta11 = convert(T, -0.45315689727309133)
    beta12 = convert(T, 0.8330937231303951)
    beta13 = convert(T, 0.3792843195533544)
    beta14 = convert(T, 0.24077885458934192)
    beta21 = convert(T, -0.4994383733810986)
    beta22 = convert(T, 0.9181786186154077)
    beta23 = convert(T, -0.25613778661003145)
    beta24 = convert(T, -0.16260245862427797)
    beta31 = convert(T, 1.4531568972730915)
    beta32 = convert(T, -0.8330937231303933)
    beta33 = convert(T, -0.3792843195533583)
    beta34 = convert(T, -0.24077885458934023)
    beta41 = convert(T, -0.4976090683622265)
    beta42 = convert(T, 0.9148155835648892)
    beta43 = convert(T, -1.4102107084476505)
    beta44 = convert(T, 0.9930041932449877)
    return FourStageSRIConstantCache(
        a021, a031, a032, a041, a042, a043, a121, a131, a132, a141,
        a142, a143, b021, b031, b032, b041, b042, b043,
        b121, b131, b132, b141, b142, b143, α1, α2,
        α3, α4, c02, c03, c04, c11, c12, c13, c14, beta11,
        beta12, beta13, beta14, beta21, beta22,
        beta23, beta24, beta31, beta32, beta33, beta34, beta41, beta42, beta43, beta44
    )
end

function alg_cache(
        alg::SOSRI2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SOSRI2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct FourStageSRICache{
        uType, randType, tabType, rateNoiseType, rateType, possibleRateType,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    chi1::randType
    chi2::randType
    chi3::randType
    tab::tabType
    g1::rateNoiseType
    g2::rateNoiseType
    g3::rateNoiseType
    g4::rateNoiseType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    E₁::rateType
    E₂::rateType
    tmp::rateType
    H02::possibleRateType
    H03::possibleRateType
end

function alg_cache(
        alg::SRIW2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi1 = copy(ΔW)
        chi2 = copy(ΔW)
        chi3 = copy(ΔW)
    else
        chi1 = zero(ΔW)
        chi2 = zero(ΔW)
        chi3 = zero(ΔW)
    end
    tab = SRIW2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    g4 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = zero(rate_prototype)
    return FourStageSRICache(
        u, uprev, chi1, chi2, chi3, tab, g1, g2, g3,
        g4, k1, k2, k3, k4, E₁, E₂, tmp, tmp, tmp
    )
end

function alg_cache(
        alg::SOSRI, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi1 = copy(ΔW)
        chi2 = copy(ΔW)
        chi3 = copy(ΔW)
    else
        chi1 = zero(ΔW)
        chi2 = zero(ΔW)
        chi3 = zero(ΔW)
    end
    tab = SOSRIConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    g4 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = zero(rate_prototype)
    return FourStageSRICache(
        u, uprev, chi1, chi2, chi3, tab, g1, g2, g3,
        g4, k1, k2, k3, k4, E₁, E₂, tmp, tmp, tmp
    )
end

function alg_cache(
        alg::SOSRI2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if ΔW isa Union{SArray, Number}
        chi1 = copy(ΔW)
        chi2 = copy(ΔW)
        chi3 = copy(ΔW)
    else
        chi1 = zero(ΔW)
        chi2 = zero(ΔW)
        chi3 = zero(ΔW)
    end
    tab = SOSRI2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
    g1 = zero(noise_rate_prototype)
    g2 = zero(noise_rate_prototype)
    g3 = zero(noise_rate_prototype)
    g4 = zero(noise_rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    E₁ = zero(rate_prototype)
    E₂ = zero(rate_prototype)
    tmp = zero(rate_prototype)
    H02 = zero(rate_prototype)
    H03 = zero(rate_prototype)
    return FourStageSRICache(
        u, uprev, chi1, chi2, chi3, tab, g1, g2, g3,
        g4, k1, k2, k3, k4, E₁, E₂, tmp, H02, H03
    )
end
