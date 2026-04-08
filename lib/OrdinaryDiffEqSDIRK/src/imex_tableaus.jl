struct ESDIRKIMEXTableau{T, T2}
    Ai::Matrix{T}
    bi::Vector{T}
    Ae::Matrix{T}
    be::Vector{T}
    c::Vector{T2}
    btilde::Union{Vector{T}, Nothing}
    ebtilde::Union{Vector{T}, Nothing}
    α::Union{Matrix{T2}, Nothing}
    order::Int
    s::Int
end

# Dispatch: each algorithm type maps to its tableau constructor
ESDIRKIMEXTableau(::ARS222, T, T2) = ARS222ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::ARS232, T, T2) = ARS232ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::ARS443, T, T2) = ARS443ESDIRKIMEXTableau(T, T2)
ESDIRKIMEXTableau(::BHR553, T, T2) = BHR553ESDIRKIMEXTableau(T, T2)

#
# ARS(2,2,2) Tableau — Ascher, Ruuth & Spiteri (1997)
# 3-stage, 2nd order, ESDIRK, ISA+GSA
#
function ARS222ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 1 - sqrt(T2(2)) / 2)
    δ = 1 - 1 / (2γ)

    s = 3

    # Implicit tableau
    Ai = zeros(T, s, s)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 2] = convert(T, 1 - γ)
    Ai[3, 3] = convert(T, γ)

    bi = zeros(T, s)
    bi[2] = convert(T, 1 - γ)
    bi[3] = convert(T, γ)

    # Explicit tableau
    Ae = zeros(T, s, s)
    Ae[2, 1] = convert(T, γ)
    Ae[3, 1] = convert(T, δ)
    Ae[3, 2] = convert(T, 1 - δ)

    be = zeros(T, s)
    be[1] = convert(T, δ)
    be[2] = convert(T, 1 - δ)

    c = zeros(T2, s)
    c[1] = zero(T2)
    c[2] = convert(T2, γ)
    c[3] = one(T2)

    return ESDIRKIMEXTableau(Ai, bi, Ae, be, c, nothing, nothing, nothing, 2, s)
end

#
# ARS(2,3,2) Tableau — Ascher, Ruuth & Spiteri (1997)
# 3-stage, 2nd order, ESDIRK, ISA (same implicit as ARS222, different explicit)
#
function ARS232ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 1 - sqrt(T2(2)) / 2)
    δ = convert(T, -2 * sqrt(T(2)) / 3)

    s = 3

    # Implicit tableau (same as ARS222)
    Ai = zeros(T, s, s)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 2] = convert(T, 1 - γ)
    Ai[3, 3] = convert(T, γ)

    bi = zeros(T, s)
    bi[2] = convert(T, 1 - γ)
    bi[3] = convert(T, γ)

    # Explicit tableau
    Ae = zeros(T, s, s)
    Ae[2, 1] = convert(T, γ)
    Ae[3, 1] = δ
    Ae[3, 2] = 1 - δ

    be = zeros(T, s)
    be[2] = convert(T, 1 - γ)
    be[3] = convert(T, γ)

    c = zeros(T2, s)
    c[1] = zero(T2)
    c[2] = convert(T2, γ)
    c[3] = one(T2)

    return ESDIRKIMEXTableau(Ai, bi, Ae, be, c, nothing, nothing, nothing, 2, s)
end

#
# ARS(4,4,3) Tableau — Ascher, Ruuth & Spiteri (1997), Table IV
# 5-stage, 3rd order, ESDIRK, ISA, γ = 1/2
#
function ARS443ESDIRKIMEXTableau(T, T2)
    γ = convert(T, 1 // 2)

    s = 5

    # Implicit tableau
    Ai = zeros(T, s, s)
    Ai[2, 2] = γ
    Ai[3, 2] = convert(T, 1 // 6)
    Ai[3, 3] = γ
    Ai[4, 2] = convert(T, -1 // 2)
    Ai[4, 3] = convert(T, 1 // 2)
    Ai[4, 4] = γ
    Ai[5, 2] = convert(T, 3 // 2)
    Ai[5, 3] = convert(T, -3 // 2)
    Ai[5, 4] = convert(T, 1 // 2)
    Ai[5, 5] = γ

    bi = zeros(T, s)
    bi[2] = convert(T, 3 // 2)
    bi[3] = convert(T, -3 // 2)
    bi[4] = convert(T, 1 // 2)
    bi[5] = convert(T, 1 // 2)

    # Explicit tableau
    Ae = zeros(T, s, s)
    Ae[2, 1] = convert(T, 1 // 2)
    Ae[3, 1] = convert(T, 11 // 18)
    Ae[3, 2] = convert(T, 1 // 18)
    Ae[4, 1] = convert(T, 5 // 6)
    Ae[4, 2] = convert(T, -5 // 6)
    Ae[4, 3] = convert(T, 1 // 2)
    Ae[5, 1] = convert(T, 1 // 4)
    Ae[5, 2] = convert(T, 7 // 4)
    Ae[5, 3] = convert(T, 3 // 4)
    Ae[5, 4] = convert(T, -7 // 4)

    # be = bi for this method
    be = zeros(T, s)
    be[2] = convert(T, 3 // 2)
    be[3] = convert(T, -3 // 2)
    be[4] = convert(T, 1 // 2)
    be[5] = convert(T, 1 // 2)

    c = zeros(T2, s)
    c[1] = zero(T2)
    c[2] = convert(T2, 1 // 2)
    c[3] = convert(T2, 2 // 3)
    c[4] = convert(T2, 1 // 2)
    c[5] = one(T2)

    return ESDIRKIMEXTableau(Ai, bi, Ae, be, c, nothing, nothing, nothing, 3, s)
end

#
# BHR(5,5,3)* Tableau — Boscarino & Russo (2009)
# 5-stage, 3rd order, ESDIRK, ISA, L-stable
#
function BHR553ESDIRKIMEXTableau(T, T2)
    γ = convert(T2, 0.435866521508460)
    b3 = 0.362863385578740
    b4 = -0.168124349878957
    c4val = 1.5
    ã53 = 1.195970114894582
    ã54 = -0.150831109536248
    a41val = 3 * c4val / 2 - c4val^2 / (4 * γ) - γ
    a43val = c4val^2 / (4 * γ) - c4val / 2
    a51val = 1 - b3 - b4 - γ
    ea41val = c4val - c4val^2 / (4 * γ)
    ea43val = c4val^2 / (4 * γ)
    ea51val = 1 + b3 - ã53 - ã54

    s = 5

    # Implicit tableau
    Ai = zeros(T, s, s)
    Ai[2, 1] = convert(T, γ)
    Ai[2, 2] = convert(T, γ)
    Ai[3, 1] = convert(T, γ)
    Ai[3, 3] = convert(T, γ)
    Ai[4, 1] = convert(T, a41val)
    Ai[4, 3] = convert(T, a43val)
    Ai[4, 4] = convert(T, γ)
    Ai[5, 1] = convert(T, a51val)
    Ai[5, 3] = convert(T, b3)
    Ai[5, 4] = convert(T, b4)
    Ai[5, 5] = convert(T, γ)

    bi = zeros(T, s)
    bi[1] = convert(T, a51val)
    bi[3] = convert(T, b3)
    bi[4] = convert(T, b4)
    bi[5] = convert(T, γ)

    # Explicit tableau
    Ae = zeros(T, s, s)
    Ae[2, 1] = convert(T, 2 * γ)
    Ae[3, 1] = convert(T, γ)
    Ae[3, 2] = convert(T, γ)
    Ae[4, 1] = convert(T, ea41val)
    Ae[4, 3] = convert(T, ea43val)
    Ae[5, 1] = convert(T, ea51val)
    Ae[5, 2] = convert(T, -b3)
    Ae[5, 3] = convert(T, ã53)
    Ae[5, 4] = convert(T, ã54)

    # be = bi for this method
    be = zeros(T, s)
    be[1] = convert(T, a51val)
    be[3] = convert(T, b3)
    be[4] = convert(T, b4)
    be[5] = convert(T, γ)

    c = zeros(T2, s)
    c[1] = zero(T2)
    c[2] = convert(T2, 2 * γ)
    c[3] = convert(T2, 2 * γ)
    c[4] = convert(T2, c4val)
    c[5] = one(T2)

    return ESDIRKIMEXTableau(Ai, bi, Ae, be, c, nothing, nothing, nothing, 3, s)
end
