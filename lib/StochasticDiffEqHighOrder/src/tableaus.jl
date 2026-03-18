"""
    RosslerSRI

Holds the Butcher tableaus for a Roessler SRI method.
"""
struct RosslerSRI{T, T2} <: Tableau
    c₀::Vector{T2}
    c₁::Vector{T2}
    A₀::Matrix{T}
    A₁::Matrix{T}
    B₀::Matrix{T}
    B₁::Matrix{T}
    α::Vector{T}
    β₁::Vector{T}
    β₂::Vector{T}
    β₃::Vector{T}
    β₄::Vector{T}
    order::Rational{Int}
end

"""
    RosslerSRA

Holds the Butcher tableaus for a Rosser SRA method.
"""
struct RosslerSRA{T, T2} <: Tableau
    c₀::Vector{T2}
    c₁::Vector{T2}
    A₀::Matrix{T}
    B₀::Matrix{T}
    α::Vector{T}
    β₁::Vector{T}
    β₂::Vector{T}
    order::Rational{Int}
end

"""
    constructSRIW1()

Constructs the tableau type for the SRIW1 method.
"""
function constructSRIW1(T = Float64, T2 = Float64)
    c₀ = [0; 3 // 4; 0; 0]
    c₁ = [0; 1 // 4; 1; 1 // 4]
    A₀ = [
        0 0 0 0
        3 // 4 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    A₁ = [
        0 0 0 0
        1 // 4 0 0 0
        1 0 0 0
        0 0 1 // 4 0
    ]
    B₀ = [
        0 0 0 0
        3 // 2 0 0 0
        0 0 0 0
        0 0 0 0
    ]
    B₁ = [
        0 0 0 0
        1 // 2 0 0 0
        -1 0 0 0
        -5 3 1 // 2 0
    ]

    α = [1 // 3; 2 // 3; 0; 0]

    β₁ = [-1; 4 // 3; 2 // 3; 0]
    β₂ = -[1; -4 // 3; 1 // 3; 0]
    β₃ = [2; -4 // 3; -2 // 3; 0]
    β₄ = [-2; 5 // 3; -2 // 3; 1]
    return RosslerSRI(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, A₁),
        map(T, B₀), map(T, B₁),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄), 3 // 2
    )
end

"""
    constructSRIW2()

Constructs the tableau type for the SRIW1 method.
"""
function constructSRIW2(T = Float64, T2 = Float64)
    c₀ = [0; 1; 1 // 2; 0]
    c₁ = [0; 1 // 4; 1; 1 // 4]
    A₀ = [
        0 0 0 0
        1 0 0 0
        1 // 4 1 // 4 0 0
        0 0 0 0
    ]
    A₁ = [
        0 0 0 0
        1 // 4 0 0 0
        1 0 0 0
        0 0 1 // 4 0
    ]
    B₀ = [
        0 0 0 0
        0 0 0 0
        1 1 // 2 0 0
        0 0 0 0
    ]
    B₁ = [
        0 0 0 0
        -1 // 2 0 0 0
        1 0 0 0
        2 -1 1 // 2 0
    ]

    α = [1 // 6; 1 // 6; 2 // 3; 0]

    β₁ = [-1; 4 // 3; 2 // 3; 0]
    β₂ = [1; -4 // 3; 1 // 3; 0]
    β₃ = [2; -4 // 3; -2 // 3; 0]
    β₄ = [-2; 5 // 3; -2 // 3; 1]
    return RosslerSRI(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, A₁),
        map(T, B₀), map(T, B₁),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄), 3 // 2
    )
end

"""
    constructSRIOpt1()

Opti6-12-11-10-01-47
"""
function constructSRIOpt1(T = Float64, T2 = Float64)
    A₀ = [
        0.0 0.0 0.0 0.0; -0.04199224421316468 0.0 0.0 0.0;
        2.842612915017106 -2.0527723684000727 0.0 0.0;
        4.338237071435815 -2.8895936137439793 2.3017575594644466 0.0
    ]
    A₁ = [
        0.0 0.0 0.0 0.0; 0.26204282091330466 0.0 0.0 0.0;
        0.20903646383505375 -0.1502377115150361 0.0 0.0;
        0.05836595312746999 0.6149440396332373 0.08535117634046772 0.0
    ]
    B₀ = [
        0.0 0.0 0.0 0.0; -0.21641093549612528 0.0 0.0 0.0;
        1.5336352863679572 0.26066223492647056 0.0 0.0;
        -1.0536037558179159 1.7015284721089472 -0.20725685784180017 0.0
    ]
    B₁ = [
        0.0 0.0 0.0 0.0; -0.5119011827621657 0.0 0.0 0.0;
        2.67767339866713 -4.9395031322250995 0.0 0.0;
        0.15580956238299215 3.2361551006624674 -1.4223118283355949 0.0
    ]
    α = [1.140099274172029, -0.6401334255743456, 0.4736296532772559, 0.026404498125060714]
    β₁ = [-1.8453464565104432, 2.688764531100726, -0.2523866501071323, 0.40896857551684956]
    β₂ = [0.4969658141589478, -0.5771202869753592, -0.12919702470322217, 0.2093514975196336]
    β₃ = [2.8453464565104425, -2.688764531100725, 0.2523866501071322, -0.40896857551684945]
    β₄ = [0.11522663875443433, -0.57877086147738, 0.2857851028163886, 0.17775911990655704]
    e = ones(size(α))
    c₀ = A₀ * e
    c₁ = A₁ * e
    return RosslerSRI(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, A₁),
        map(T, B₀), map(T, B₁),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄), 3 // 2
    )
end

"""
    constructSRIOpt2()

Opti6-12-11-10-01-47
"""
function constructSRIOpt2(T = Float64, T2 = Float64)
    c0 = [0.0, 0.13804532298278663, 0.9999999999999992, 0.9999999999999994]
    c1 = [0.0, 0.45605532163856893, 0.999999999999996, 0.9999999999999962]
    A0 = [
        0.0 0.0 0.0 0.0
        0.13804532298278663 0.0 0.0 0.0
        0.5818361298250374 0.4181638701749618 0.0 0.0
        0.4670018408674211 0.8046204792187386 -0.27162232008616016 0.0
    ]
    A1 = [
        0.0 0.0 0.0 0.0
        0.45605532163856893 0.0 0.0 0.0
        0.7555807846451692 0.24441921535482677 0.0 0.0
        0.6981181143266059 0.3453277086024727 -0.04344582292908241 0.0
    ]
    B0 = [
        0.0 0.0 0.0 0.0
        0.08852381537667678 0.0 0.0 0.0
        1.0317752458971061 0.4563552922077882 0.0 0.0
        1.73078280444124 -0.46089678470929774 -0.9637509618944188 0.0
    ]
    B1 = [
        0.0 0.0 0.0 0.0
        0.6753186815412179 0.0 0.0 0.0
        -0.07452812525785148 -0.49783736486149366 0.0 0.0
        -0.5591906709928903 0.022696571806569924 -0.8984927888368557 0.0
    ]
    α = [-0.15036858140642623, 0.7545275856696072, 0.686995463807979, -0.2911544680711602]
    β1 = [-0.45315689727309133, 0.8330937231303951, 0.3792843195533544, 0.24077885458934192]
    β2 = [
        -0.4994383733810986, 0.9181786186154077, -0.25613778661003145, -0.16260245862427797,
    ]
    β3 = [
        1.4531568972730915, -0.8330937231303933, -0.3792843195533583, -0.24077885458934023,
    ]
    β4 = [-0.4976090683622265, 0.9148155835648892, -1.4102107084476505, 0.9930041932449877]
    return RosslerSRI(
        map(T2, c0), map(T2, c1),
        map(T, A0), map(T, A1),
        map(T, B0), map(T, B1),
        map(T, α), map(T, β1), map(T, β2),
        map(T, β3), map(T, β4), 3 // 2
    )
end

"""
    constructSRA1()

Constructs the taleau type for the SRA1 method.
"""
function constructSRA1(T = Float64, T2 = Float64)
    α = [1 // 3; 2 // 3]
    β₁ = [1; 0]
    β₂ = [-1; 1]
    A₀ = [
        0 0
        3 // 4 0
    ]
    B₀ = [
        0 0
        3 // 2 0
    ]
    c₀ = [0; 3 // 4]
    c₁ = [1; 0]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end

"""
    constructSRA2()

Constructs the taleau type for the SRA2 method.
"""
function constructSRA2(T = Float64, T2 = Float64)
    α = [1 // 3; 2 // 3]
    β₁ = [0; 1]
    β₂ = [3 // 2; -3 // 2]
    A₀ = [
        0 0
        3 // 4 0
    ]
    B₀ = [
        0 0
        3 // 2 0
    ]
    c₀ = [0; 3 // 4]
    c₁ = [1 // 3; 1]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end

"""
    constructSRA3()

Constructs the taleau type for the SRA3 method.
"""
function constructSRA3(T = Float64, T2 = Float64)
    α = [1 // 6; 1 // 6; 2 // 3]
    β₁ = [1; 0; 0]
    β₂ = [-1; 1; 0]
    A₀ = [
        0 0 0
        1 0 0
        1 // 4 1 // 4 0
    ]
    B₀ = [
        0 0 0
        0 0 0
        1 1 // 2 0
    ]
    c₀ = [0; 1; 1 // 2]
    c₁ = [1; 0; 0]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end

"""
    constructSOSRA()

Constructs the taleau type for the SOSRA method.
"""
function constructSOSRA(T = Float64, T2 = Float64)
    a1 = 0.2889874966892885
    a2 = 0.6859880440839937
    a3 = 0.025024459226717772
    c01 = 0
    c02 = 0.6923962376159507
    c03 = 1
    c11 = 0
    c12 = 0.041248171110700504
    c13 = 1
    b11 = -16.792534242221663
    b12 = 17.514995785380226
    b13 = 0.27753845684143835
    b21 = 0.4237535769069274
    b22 = 0.6010381474428539
    b23 = -1.0247917243497813
    A021 = 0.6923962376159507
    A031 = -3.1609142252828395
    A032 = 4.1609142252828395
    B021 = 1.3371632704399763
    B031 = 1.442371048468624
    B032 = 1.8632741501139225
    α = [a1; a2; a3]
    β₁ = [b11; b12; b13]
    β₂ = [b21; b22; b23]
    A₀ = [
        0 0 0
        A021 0 0
        A031 A032 0
    ]
    B₀ = [
        0 0 0
        B021 0 0
        B031 B032 0
    ]
    c₀ = [c01; c02; c03]
    c₁ = [c11; c12; c13]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end

"""
    constructSOSRA2()

Constructs the taleau type for the SOSRA method.
"""
function constructSOSRA2(T = Float64, T2 = Float64)
    a1 = 0.4999999999999998
    a2 = -0.9683897375354181
    a3 = 1.4683897375354185
    c01 = 0
    c02 = 1
    c03 = 1
    c11 = 0
    c12 = 1.0
    c13 = 1
    b11 = 0.0
    b12 = 0.92438032145683
    b13 = 0.07561967854316998
    b21 = 1.0
    b22 = -0.8169981105823436
    b23 = -0.18300188941765633
    A021 = 1
    A031 = 0.9511849235504364
    A032 = 0.04881507644956362
    B021 = 0.7686101171003622
    B031 = 0.43886792994934987
    B032 = 0.7490415909204886
    α = [a1; a2; a3]
    β₁ = [b11; b12; b13]
    β₂ = [b21; b22; b23]
    A₀ = [
        0 0 0
        A021 0 0
        A031 A032 0
    ]
    B₀ = [
        0 0 0
        B021 0 0
        B031 B032 0
    ]
    c₀ = [c01; c02; c03]
    c₁ = [c11; c12; c13]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end

"""
    checkSRIOrder(RosslerSRI)

Determines whether the order conditions are met via the tableaus of the SRI method.
"""
function checkSRIOrder(RosslerSRI; tol = 1.0e-6)
    (; c₀, c₁, A₀, A₁, B₀, B₁, α, β₁, β₂, β₃, β₄) = RosslerSRI
    e = ones(size(α))
    conditions = Vector{Bool}(undef, 25)
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₁, e) - 1) < tol
    conditions[3] = abs(dot(β₂, e) - 0) < tol
    conditions[4] = abs(dot(β₃, e) - 0) < tol
    conditions[5] = abs(dot(β₄, e) - 0) < tol
    conditions[6] = abs(sum(β₁' * B₁) - 0) < tol
    conditions[7] = abs(sum(β₂' * B₁) - 1) < tol
    conditions[8] = abs(sum(β₃' * B₁) - 0) < tol
    conditions[9] = abs(sum(β₄' * B₁) - 0) < tol
    conditions[10] = abs(sum(α' * A₀) - 0.5) < tol
    conditions[11] = abs(sum(α' * B₀) - 1) < tol
    conditions[12] = abs(sum(α' * (B₀ * e) .^ 2) - 3 / 2) < tol
    conditions[13] = abs(sum(β₁' * A₁) - 1) < tol
    conditions[14] = abs(sum(β₂' * A₁) - 0) < tol
    conditions[15] = abs(sum(β₃' * A₁) + 1) < tol
    conditions[16] = abs(sum(β₄' * A₁) - 0) < tol
    conditions[17] = abs(sum(β₁' * (B₁ * e) .^ 2) - 1) < tol
    conditions[18] = abs(sum(β₂' * (B₁ * e) .^ 2) - 0) < tol
    conditions[19] = abs(sum(β₃' * (B₁ * e) .^ 2) + 1) < tol
    conditions[20] = abs(sum(β₄' * (B₁ * e) .^ 2) - 2) < tol
    conditions[22] = abs(sum(β₂' * B₁ * (B₁ * e)) - 0) < tol
    conditions[21] = abs(sum(β₁' * B₁ * (B₁ * e)) - 0) < tol
    conditions[23] = abs(sum(β₃' * B₁ * (B₁ * e)) - 0) < tol
    conditions[24] = abs(sum(β₄' * B₁ * (B₁ * e)) - 1) < tol
    conditions[25] = abs.(0.5 * β₁' * (A₁ * (B₀ * e)) + (1 / 3) * β₃' * (A₁ * (B₀ * e)))[1] .< tol
    return (conditions)
end

"""
    checkSRAOrder(RosslerSRI)

Determines whether the order conditions are met via the tableaus of the SRA method.
"""
function checkSRAOrder(SRA; tol = 1.0e-6)
    (; c₀, c₁, A₀, B₀, α, β₁, β₂) = SRA
    e = ones(size(α))
    conditions = Vector{Bool}(undef, 8)
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₁, e) - 1) < tol
    conditions[3] = abs(dot(β₂, e) - 0) < tol
    conditions[4] = abs(dot(α, B₀ * e) - 1) .< tol
    conditions[5] = abs(dot(α, A₀ * e) - 1 / 2) .< tol
    conditions[6] = abs(dot(α, (B₀ * e) .^ 2) - 3 / 2) .< tol
    conditions[7] = abs(dot(β₁, c₁) - 1) < tol
    conditions[8] = abs(dot(β₂, c₁) + 1) < tol
    return (conditions)
end

"""
    constructSKenCarp()

Constructs the tableau type for the implicit SKenCarp method as a RosslerSRA tableau.
"""
function constructSKenCarp(T = Float64, T2 = Float64)
    γ = convert(T, 0.435866521508459)
    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)
    # bhat1 = convert(T,2756255671327//12835298489170)
    # bhat2 = -convert(T,10771552573575//22201958757719)
    # bhat3 = convert(T,9247589265047//10645013368117)
    # bhat4 = convert(T,2193209047091//5459859503100)
    btilde1 = convert(T, 0.027099261876665316) # bhat1-a41
    btilde2 = convert(T, 0.11013520969201586) # bhat2-a42
    btilde3 = convert(T, -0.10306492520138458) # bhat3-a43
    btilde4 = convert(T, -0.0341695463672966) # bhat4-γ
    c3 = convert(T2, 0.6)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    nb021 = convert(T, -12.246764387585056)
    nb043 = convert(T, -14.432096958608753)
    α = [a41; a42; a43; γ]
    β₁ = [0, 0, 0, 1]
    β₂ = [1, 0, 0, -1]
    A₀ = [
        0 0 0 0
        γ γ 0 0
        a31 a32 γ 0
        a41 a42 a43 γ
    ]
    B₀ = [
        0 0 0 0
        nb021 0 0 0
        0 0 0 0
        0 0 nb043 0
    ]
    c₀ = [0, c2, c3, 1]
    c₁ = [0, 0, 0, 1]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end

"""
    constructExplicitSKenCarp()

Constructs the tableau type for the explicit part of SKenCarp as a RosslerSRA tableau.
"""
function constructExplicitSKenCarp(T = Float64, T2 = Float64)
    γ = convert(T, 0.435866521508459)
    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)
    # bhat1 = convert(T,2756255671327//12835298489170)
    # bhat2 = -convert(T,10771552573575//22201958757719)
    # bhat3 = convert(T,9247589265047//10645013368117)
    # bhat4 = convert(T,2193209047091//5459859503100)
    btilde1 = convert(T, 0.027099261876665316) # bhat1-a41
    btilde2 = convert(T, 0.11013520969201586) # bhat2-a42
    btilde3 = convert(T, -0.10306492520138458) # bhat3-a43
    btilde4 = convert(T, -0.0341695463672966) # bhat4-γ
    c3 = convert(T2, 0.6)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    ea21 = convert(T, 1767732205903 // 2027836641118)
    ea31 = convert(T, 5535828885825 // 10492691773637)
    ea32 = convert(T, 788022342437 // 10882634858940)
    ea41 = convert(T, 6485989280629 // 16251701735622)
    ea42 = -convert(T, 4246266847089 // 9704473918619)
    ea43 = convert(T, 10755448449292 // 10357097424841)
    eb1 = convert(T, 1471266399579 // 7840856788654)
    eb2 = convert(T, -4482444167858 // 7529755066697)
    eb3 = convert(T, 11266239266428 // 11593286722821)
    eb4 = convert(T, 1767732205903 // 4055673282236)

    nb021 = convert(T, -12.246764387585056)
    nb043 = convert(T, -14.432096958608753)

    α = [eb1; eb2; eb3; eb4]
    β₁ = [0, 0, 0, 1]
    β₂ = [1, 0, 0, -1]
    A₀ = [
        0 0 0 0
        ea21 0 0 0
        ea31 ea32 0 0
        ea41 ea42 ea43 0
    ]
    B₀ = [
        0 0 0 0
        nb021 0 0 0
        0 0 0 0
        0 0 nb043 0
    ]
    c₀ = [0, c2, c3, 1]
    c₁ = [0, 0, 0, 1]
    return RosslerSRA(
        map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2
    )
end
