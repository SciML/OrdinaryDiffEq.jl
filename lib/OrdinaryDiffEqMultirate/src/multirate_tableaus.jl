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

# ── MRI-GARK: explicit infinitesimal GARK ─────────────────────────────────────

struct MRIGARKExplicitTableau{T}
    Γ::Vector{T}  # diagonal fast rates, length s
    γ::Matrix{T}  # lower-triangular slow couplings, s×s
    c::Vector{T}  # slow abscissae for fS evaluations, length s
end

function MRIGARKERK22aTableau(::Type{T}) where {T}
    Γ = T[1 // 2, 1 // 2]
    γ = T[1 // 2 0; -1 // 2 1]
    c = T[0, 1]
    return MRIGARKExplicitTableau{T}(Γ, γ, c)
end

function MRIGARKERK22bTableau(::Type{T}) where {T}
    Γ = T[1, 0]
    γ = T[1 0; -1 // 2 1 // 2]
    c = T[0, 1]
    return MRIGARKExplicitTableau{T}(Γ, γ, c)
end

mri_gark_tableau(::MRIGARKERK22a, ::Type{T}) where {T} = MRIGARKERK22aTableau(T)
mri_gark_tableau(::MRIGARKERK22b, ::Type{T}) where {T} = MRIGARKERK22bTableau(T)

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
