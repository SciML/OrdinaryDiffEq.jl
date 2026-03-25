struct MronRoe <: AbstractIteratedIntegralAlgorithm end

convorder(::MronRoe) = 1 // 1
errcoeff(m, h, ::MronRoe, ::MaxLp{2}) = √m * h / (√12 * π)
norv(m, n, ::MronRoe) = 2 * m * n + (m^2 + m) ÷ 2

function levyarea(W::AbstractVector{T}, n::Integer, alg::MronRoe;
        rng::AbstractRNG = default_rng()) where {T <: AbstractFloat}
    m = length(W)
    # MronRoe uses transposed layout: X is m×n, Y is n×m
    X = randn(rng, T, m, n)
    Y = randn(rng, T, n, m)
    Ψ = randn(rng, T, m)
    n_tail = (m^2 - m) ÷ 2
    tail_extra = randn(rng, T, n_tail)
    return _mronroe_area(W, X, Y, Ψ, tail_extra, n, m)
end

function levyarea(W::AbstractVector{T}, n::Integer, alg::MronRoe,
        coeffs::LevyAreaCoefficients{T}) where {T <: AbstractFloat}
    m = coeffs.m
    Ψ = coeffs.tail[1:m]
    tail_extra = coeffs.tail[(m + 1):end]
    return _mronroe_area(W, coeffs.X, copy(coeffs.Y), Ψ, tail_extra, n, m)
end

function _mronroe_area(W::AbstractVector{T}, X, Y, Ψ, tail_extra, n, m) where {T}
    Y = copy(Y)
    Y .= (Y .- √(T(2)) .* W') ./ (1:n)
    A = X * Y

    a = T(sqrt(2 * trigamma(n + 1)))
    A .+= a .* W .* Ψ'

    idx = 0
    for i in 1:m
        @inbounds A[i, i] = zero(T)
        for j in (i + 1):m
            idx += 1
            @inbounds A[i, j] = (A[i, j] + a * tail_extra[idx] - A[j, i]) / (2 * T(π))
            @inbounds A[j, i] = -A[i, j]
        end
    end

    return A
end
