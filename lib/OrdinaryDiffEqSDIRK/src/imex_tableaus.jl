struct IMEXTableau{T, T2}
    Ai::Matrix{T}
    bi::Vector{T}
    Ae::Matrix{T}
    be::Vector{T}
    c::Vector{T2}
    btilde::Vector{T}
    ebtilde::Union{Vector{T}, Nothing}
    α::Union{Matrix{T2}, Nothing}
    order::Int
    s::Int
end

function KenCarp3IMEXTableau(T, T2)
    γ = convert(T, 1767732205903 // 4055673282236)

    a31 = convert(T, 2746238789719 // 10658868560708)
    a32 = -convert(T, 640167445237 // 6845629431997)
    a41 = convert(T, 1471266399579 // 7840856788654)
    a42 = -convert(T, 4482444167858 // 7529755066697)
    a43 = convert(T, 11266239266428 // 11593286722821)

    btilde1 = convert(
        T,
        BigInt(681815649026867975666107) //
            BigInt(25159934323302256049469295)
    )
    btilde2 = convert(
        T,
        BigInt(18411887981491912264464127) //
            BigInt(167175311446532472108584143)
    )
    btilde3 = -convert(
        T,
        BigInt(12719313754959329011138489) //
            BigInt(123410692144842870217698057)
    )
    btilde4 = -convert(
        T,
        BigInt(47289384293135913063989) //
            BigInt(1383962894467812063558225)
    )

    c3 = convert(T2, 3 // 5)
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
    ebtilde1 = convert(
        T,
        BigInt(681815649026867975666107) //
            BigInt(25159934323302256049469295)
    )
    ebtilde2 = convert(
        T,
        BigInt(18411887981491912264464127) //
            BigInt(167175311446532472108584143)
    )
    ebtilde3 = -convert(
        T,
        BigInt(12719313754959329011138489) //
            BigInt(123410692144842870217698057)
    )
    ebtilde4 = -convert(
        T,
        BigInt(47289384293135913063989) //
            BigInt(1383962894467812063558225)
    )

    s = 4
    Ai = zeros(T, s, s)
    Ai[2, 1] = γ
    Ai[2, 2] = γ
    Ai[3, 1] = a31
    Ai[3, 2] = a32
    Ai[3, 3] = γ
    Ai[4, 1] = a41
    Ai[4, 2] = a42
    Ai[4, 3] = a43
    Ai[4, 4] = γ

    bi_vec = zeros(T, s)
    bi_vec[1] = a41
    bi_vec[2] = a42
    bi_vec[3] = a43
    bi_vec[4] = γ

    Ae = zeros(T, s, s)
    Ae[2, 1] = ea21
    Ae[3, 1] = ea31
    Ae[3, 2] = ea32
    Ae[4, 1] = ea41
    Ae[4, 2] = ea42
    Ae[4, 3] = ea43

    be_vec = zeros(T, s)
    be_vec[1] = eb1
    be_vec[2] = eb2
    be_vec[3] = eb3
    be_vec[4] = eb4

    c_vec = zeros(T2, s)
    c_vec[1] = zero(T2)
    c_vec[2] = convert(T2, 2γ)
    c_vec[3] = c3
    c_vec[4] = one(T2)

    btilde_vec = zeros(T, s)
    btilde_vec[1] = btilde1
    btilde_vec[2] = btilde2
    btilde_vec[3] = btilde3
    btilde_vec[4] = btilde4

    ebtilde_vec = zeros(T, s)
    ebtilde_vec[1] = ebtilde1
    ebtilde_vec[2] = ebtilde2
    ebtilde_vec[3] = ebtilde3
    ebtilde_vec[4] = ebtilde4

    α_mat = zeros(T2, s, s)
    α_mat[3, 1] = α31
    α_mat[3, 2] = α32
    α_mat[4, 1] = α41
    α_mat[4, 2] = α42

    return IMEXTableau(Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 3, s)
end

function ARS343Tableau(T, T2)
    γ = convert(T, 4358665215084590 // 10000000000000000)

    s = 4

    c2 = γ
    c3 = (one(T2) + convert(T2, γ)) / 2
    c4 = one(T2)

    a32_i = (one(T) - γ) / 2

    b3_i = (one(T) / 2 - 2γ + γ^2) / ((one(T) - γ) / 2)
    b2_i = one(T) - γ - b3_i

    Ai = zeros(T, s, s)
    Ai[2, 2] = γ
    Ai[3, 2] = a32_i
    Ai[3, 3] = γ
    Ai[4, 2] = b2_i
    Ai[4, 3] = b3_i
    Ai[4, 4] = γ

    bi_vec = T[zero(T), b2_i, b3_i, γ]

    ae21 = γ
    ae31 = convert(T, 3212788860 // 10000000000)
    ae32 = convert(T, 3966543748 // 10000000000)
    ae41 = -convert(T, 1058582960 // 10000000000)
    ae42 = convert(T, 5529291479 // 10000000000)
    ae43 = convert(T, 5529291479 // 10000000000)

    Ae = zeros(T, s, s)
    Ae[2, 1] = ae21
    Ae[3, 1] = ae31
    Ae[3, 2] = ae32
    Ae[4, 1] = ae41
    Ae[4, 2] = ae42
    Ae[4, 3] = ae43

    be_vec = T[zero(T), b2_i, b3_i, γ]

    c_vec = T2[zero(T2), convert(T2, c2), c3, c4]

    btilde_vec = bi_vec .- T[zero(T), zero(T), zero(T), one(T)]
    ebtilde_vec = be_vec .- T[zero(T), zero(T), zero(T), one(T)]

    α_mat = zeros(T2, s, s)

    return IMEXTableau(Ai, bi_vec, Ae, be_vec, c_vec,
        btilde_vec, ebtilde_vec, α_mat, 3, s)
end
