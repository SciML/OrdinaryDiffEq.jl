struct MRIGARKExplicit22Tableau{T}
    Γ11::T  # fast rate, stage 1
    Γ22::T  # fast rate, stage 2 (0 when c₂ = 1)
    γ11::T  # f2(u_n) coupling, stage 1
    γ21::T  # f2(u_n) coupling, stage 2
    γ22::T  # f2(Y₂) coupling, stage 2
end

function MRIGARKERK22aTableau(::Type{T}) where {T}
    return MRIGARKExplicit22Tableau{T}(
        T(1 // 2), T(1 // 2), T(1 // 2), T(-1 // 2), T(1)
    )
end

function MRIGARKERK22bTableau(::Type{T}) where {T}
    return MRIGARKExplicit22Tableau{T}(
        T(1), T(0), T(1), T(-1 // 2), T(1 // 2)
    )
end
