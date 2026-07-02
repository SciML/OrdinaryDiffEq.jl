# ── MRAB: Adams–Bashforth weights ─────────────────────────────────────────────

struct MRABTableau{T}
    β::Vector{Vector{T}}  # β[j] = explicit Adams–Bashforth-j weights, β[j][1] for the newest rate
end

function MRABTableau(k::Int, ::Type{T}) where {T}
    1 <= k <= 5 ||
        throw(ArgumentError("MRAB: unsupported Adams order k=$k (supported: 1..5)"))
    allβ = [
        T[1],
        T[3 // 2, -1 // 2],
        T[23 // 12, -16 // 12, 5 // 12],
        T[55 // 24, -59 // 24, 37 // 24, -9 // 24],
        T[1901 // 720, -2774 // 720, 2616 // 720, -1274 // 720, 251 // 720],
    ]
    return MRABTableau{T}(allβ[1:k])
end

# ── MRI-GARK: explicit infinitesimal GARK (Sandu 2019) ────────────────────────
# Stage i integrates the fast IVP v' = Δcᵢ·dt·f1(v) + dt·Σⱼ ωᵢⱼ(τ)·f2(z_{j-1})
# over normalized τ∈[0,1], where ωᵢⱼ(τ) = W0[i,j] + W1[i,j]·τ. Inner solver is an
# explicit RK of order `q`. Wemb0/Wemb1 (length 0 if absent) give the embedded
# last-stage coupling ωᵉᵐᵇⱼ(τ) = Wemb0[j] + Wemb1[j]·τ for the error estimate.

struct MRIGARKTableau{T}
    Δc::Vector{T}     # fast sub-interval length per stage, length s
    W0::Matrix{T}     # constant coupling, s×s lower-triangular
    W1::Matrix{T}     # linear-in-τ coupling, s×s lower-triangular
    Wemb0::Vector{T}  # embedded last-stage constant coupling (length s) or empty
    Wemb1::Vector{T}  # embedded last-stage linear-in-τ coupling, or empty (constant)
    γ0::Vector{T}     # implicit endpoint coupling per stage (0 ⇒ explicit stage), length s
    q::Int            # inner RK order
end

function MRIGARKERK22aTableau(::Type{T}) where {T}
    Δc = T[1 // 2, 1 // 2]
    W0 = T[1 // 2 0; -1 // 2 1]
    W1 = zeros(T, 2, 2)
    return MRIGARKTableau{T}(Δc, W0, W1, T[], T[], zeros(T, 2), 2)
end

function MRIGARKERK22bTableau(::Type{T}) where {T}
    Δc = T[1, 0]
    W0 = T[1 0; -1 // 2 1 // 2]
    W1 = zeros(T, 2, 2)
    return MRIGARKTableau{T}(Δc, W0, W1, T[], T[], zeros(T, 2), 2)
end

function MRIGARKERK33aTableau(::Type{T}) where {T}
    Δc = T[1 // 3, 1 // 3, 1 // 3]
    W0 = T[1 // 3 0 0; -1 // 3 2 // 3 0; 0 -2 // 3 1]
    W1 = T[0 0 0; 0 0 0; 1 // 2 0 -1 // 2]
    Wemb0 = T[1 // 12, -1 // 3, 7 // 12]
    return MRIGARKTableau{T}(Δc, W0, W1, Wemb0, T[], zeros(T, 3), 3)
end

function MRIGARKERK45aTableau(::Type{T}) where {T}
    Δc = T[1 // 5, 1 // 5, 1 // 5, 1 // 5, 1 // 5]
    W0 = zeros(T, 5, 5)
    W0[1, 1] = 1 // 5
    W0[2, 1] = -53 // 16; W0[2, 2] = 281 // 80
    W0[3, 1] = -36562993 // 71394880; W0[3, 2] = 34903117 // 17848720
    W0[3, 3] = -88770499 // 71394880
    W0[4, 1] = -7631593 // 71394880; W0[4, 2] = -166232021 // 35697440
    W0[4, 3] = 6068517 // 1519040; W0[4, 4] = 8644289 // 8924360
    W0[5, 1] = 277061 // 303808; W0[5, 2] = -209323 // 1139280
    W0[5, 3] = -1360217 // 1139280; W0[5, 4] = -148789 // 56964; W0[5, 5] = 147889 // 45120
    W1 = zeros(T, 5, 5)
    W1[2, 1] = 503 // 80; W1[2, 2] = -503 // 80
    W1[3, 1] = -1365537 // 35697440; W1[3, 2] = 4963773 // 7139488
    W1[3, 3] = -1465833 // 2231090
    W1[4, 1] = 66974357 // 35697440; W1[4, 2] = 21445367 // 7139488
    W1[4, 3] = -3; W1[4, 4] = -8388609 // 4462180
    W1[5, 1] = -18227 // 7520; W1[5, 2] = 2; W1[5, 3] = 1; W1[5, 4] = 5
    W1[5, 5] = -41933 // 7520
    Wemb0 = T[
        -88227 // 47470, 756870829 // 340217490, -713704111 // 1360869960,
        -31967827 // 340217490, 129673 // 286680,
    ]
    Wemb1 = T[6213 // 1880, -6213 // 1880, 0, 0, 0]
    return MRIGARKTableau{T}(Δc, W0, W1, Wemb0, Wemb1, zeros(T, 5), 4)
end

function MRIGARKIRK21aTableau(::Type{T}) where {T}
    Δc = T[1, 0]
    W0 = T[1 0; -1 // 2 0]
    W1 = zeros(T, 2, 2)
    γ0 = T[0, 1 // 2]
    return MRIGARKTableau{T}(Δc, W0, W1, T[], T[], γ0, 2)
end

function MRIGARKESDIRK34aTableau(::Type{T}) where {T}
    β = T(0.4358665215084589994160194511935568425)
    Δc = T[1//3, 0, 1//3, 0, 1//3, 0]
    W0 = zeros(T, 6, 6)
    W0[1, 1] = 1 // 3
    W0[2, 1] = -β
    W0[3, 1] = -0.3045790611944504970424837655380884888
    W0[3, 3] = 0.6379123945277838303758170988714218222
    W0[4, 1] = 0.2116913105640266601676536489364004869
    W0[4, 3] = -0.6475578320724856595836731001299573294
    W0[5, 1] = 0.4454209388055495029575162344619115112
    W0[5, 3] = 0.8813784805616198280398949036456491923
    W0[5, 5] = -0.9934660860338359976640778047742273701
    W0[6, 1] = -β
    W1 = zeros(T, 6, 6)
    γ0 = T[0, β, 0, β, 0, β]
    return MRIGARKTableau{T}(Δc, W0, W1, T[], T[], γ0, 3)
end

mri_gark_tableau(::MRIGARKERK22a, ::Type{T}) where {T} = MRIGARKERK22aTableau(T)
mri_gark_tableau(::MRIGARKERK22b, ::Type{T}) where {T} = MRIGARKERK22bTableau(T)
mri_gark_tableau(::MRIGARKERK33a, ::Type{T}) where {T} = MRIGARKERK33aTableau(T)
mri_gark_tableau(::MRIGARKERK45a, ::Type{T}) where {T} = MRIGARKERK45aTableau(T)
mri_gark_tableau(::MRIGARKIRK21a, ::Type{T}) where {T} = MRIGARKIRK21aTableau(T)
mri_gark_tableau(::MRIGARKESDIRK34a, ::Type{T}) where {T} = MRIGARKESDIRK34aTableau(T)

# ── MIS: multirate infinitesimal step ─────────────────────────────────────────

struct MISTableau{T}
    α::Matrix{T}   # strictly lower-triangular stage couplings
    β::Matrix{T}   # strictly lower-triangular slow-rate couplings
    γ::Matrix{T}   # strictly lower-triangular correction couplings
    d::Vector{T}   # inner-interval lengths, d_i = Σ_j β[i,j]
    c::Vector{T}   # stage abscissae, (I − α − γ)⁻¹ d
    ctilde::Vector{T}  # inner start abscissae, α c
end

function MIS2Tableau(::Type{T}) where {T}
    α = T[
        0 0 0 0
        0 0 0 0
        0 0.53694656671 0 0
        0 0.480892968551 0.500561163566 0
    ]
    β = T[
        0 0 0 0
        0.126848494553 0 0 0
        -0.784838278826 1.37442675268 0 0
        -0.0456727081749 -0.0087508227119 0.524775788629 0
    ]
    γ = T[
        0 0 0 0
        0 0 0 0
        0 0.652465126004 0 0
        0 -0.0732769849457 0.14490243042 0
    ]
    d = T[sum(@view β[i, :]) for i in axes(β, 1)]
    c = (LinearAlgebra.I(length(d)) - α - γ) \ d
    ctilde = α * c
    return MISTableau{T}(α, β, γ, d, c, ctilde)
end
