struct MRIGARKExplicit22Tableau{T}
    Γ11::T
    Γ22::T  # 0 for ERK22b (c₂=1); no inner fast loop in stage 2
    γ11::T
    γ21::T
    γ22::T
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
