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
