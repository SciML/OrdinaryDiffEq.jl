struct Fourier <: AbstractIteratedIntegralAlgorithm end

convorder(::Fourier) = 1 // 2
errcoeff(m, h, ::Fourier, ::MaxLp{2}) = h * √3 / (√2 * π)
norv(m, n, ::Fourier) = 2 * m * n

function levyarea(W::AbstractVector{T}, n::Integer, alg::Fourier;
        rng::AbstractRNG = default_rng()) where {T <: AbstractFloat}
    m = length(W)
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    return _fourier_area(W, X, Y, n, m)
end

function levyarea(W::AbstractVector{T}, n::Integer, alg::Fourier,
        coeffs::LevyAreaCoefficients{T}) where {T <: AbstractFloat}
    return _fourier_area(W, coeffs.X, copy(coeffs.Y), n, coeffs.m)
end

function _fourier_area(W::AbstractVector{T}, X, Y, n, m) where {T}
    Y = copy(Y)
    Y .= (Y .- √(T(2)) .* W) ./ (1:n)'
    A = Y * X
    for i in 1:m
        @inbounds A[i, i] = zero(T)
        for j in (i + 1):m
            @inbounds A[i, j] = (A[i, j] - A[j, i]) / (2 * T(π))
            @inbounds A[j, i] = -A[i, j]
        end
    end
    return A
end
