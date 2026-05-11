const CompiledFloats = Union{Float32, Float64}

struct SKenCarpTableau{T, T2}
    γ::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    c3::T2
    α31::T
    α32::T
    α41::T
    α42::T
    ea21::T
    ea31::T
    ea32::T
    ea41::T
    ea42::T
    ea43::T
    eb1::T
    eb2::T
    eb3::T
    eb4::T
    ebtilde1::T
    ebtilde2::T
    ebtilde3::T
    ebtilde4::T
    nb021::T
    nb043::T
end

#=
# KenCarp3
# Predict z4 from Hermite z2 and z1
# Not z3 because c3 < c2 !

θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
θ = c4/c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
=#
function SKenCarpTableau(::Type{T}, ::Type{T2}) where {
        T <: CompiledFloats, T2 <: CompiledFloats,
    }
    γ = convert(T, 0.435866521508459)
    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)
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

    ea21 = convert(T, 0.871733043016918)
    ea31 = convert(T, 0.5275890119763004)
    ea32 = convert(T, 0.0724109880236996)
    ea41 = convert(T, 0.3990960076760701)
    ea42 = -convert(T, 0.4375576546135194)
    ea43 = convert(T, 1.0384616469374492)
    eb1 = convert(T, 0.18764102434672383)
    eb2 = convert(T, -0.595297473576955)
    eb3 = convert(T, 0.9717899277217721)
    eb4 = convert(T, 0.435866521508459)
    ebtilde1 = convert(T, 0.027099261876665316)
    ebtilde2 = convert(T, 0.11013520969201586)
    ebtilde3 = -convert(T, 0.10306492520138458)
    ebtilde4 = -convert(T, 0.0341695463672966)

    nb021 = convert(T, -12.246764387585056)
    nb043 = convert(T, -14.432096958608753)
    return SKenCarpTableau(
        γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31,
        α32, α41, α42, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4,
        ebtilde1, ebtilde2, ebtilde3, ebtilde4, nb021, nb043
    )
end

function SKenCarpTableau(::Type{T}, ::Type{T2}) where {T, T2}
    γ = convert(T, 1767732205903 // 4055673282236)
    a31 = convert(T, 2746238789719 // 10658868560708)
    a32 = -convert(T, 640167445237 // 6845629431997)
    a41 = convert(T, 1471266399579 // 7840856788654)
    a42 = -convert(T, 4482444167858 // 7529755066697)
    a43 = convert(T, 11266239266428 // 11593286722821)
    btilde1 = convert(T, BigInt(681815649026867975666107) // BigInt(25159934323302256049469295)) # bhat1-a41
    btilde2 = convert(T, BigInt(18411887981491912264464127) // BigInt(167175311446532472108584143)) # bhat2-a42
    btilde3 = convert(T, BigInt(-12719313754959329011138489) // BigInt(123410692144842870217698057)) # bhat3-a43
    btilde4 = convert(T, BigInt(-47289384293135913063989) // BigInt(1383962894467812063558225)) # bhat4-γ
    c3 = convert(T2, 3 // 5)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    # Explicit Tableau
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
    ebtilde1 = convert(T, BigInt(681815649026867975666107) // BigInt(25159934323302256049469295))
    ebtilde2 = convert(T, BigInt(18411887981491912264464127) // BigInt(167175311446532472108584143))
    ebtilde3 = -convert(T, BigInt(12719313754959329011138489) // BigInt(123410692144842870217698057))
    ebtilde4 = -convert(T, BigInt(47289384293135913063989) // BigInt(1383962894467812063558225))

    # Noise Tableau
    nb021 = convert(
        T,
        parse(
            BigFloat,
            "-12.246764387585055918338744103409192607986567514699471403397969732723452087723101"
        )
    )
    nb043 = convert(
        T,
        parse(
            BigFloat,
            "-14.432096958608752822047165680776748797565142459789556194474191884258734697161106"
        )
    )
    return SKenCarpTableau(
        γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31,
        α32, α41, α42, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4,
        ebtilde1, ebtilde2, ebtilde3, ebtilde4, nb021, nb043
    )
end
