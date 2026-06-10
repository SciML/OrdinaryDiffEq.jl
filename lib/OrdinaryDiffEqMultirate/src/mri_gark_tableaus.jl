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
