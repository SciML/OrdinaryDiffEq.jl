struct Milstein <: AbstractIteratedIntegralAlgorithm end

convorder(::Milstein) = 1 // 2
errcoeff(m, h, ::Milstein, ::MaxLp{2}) = h / (√2 * π)
norv(m, n, ::Milstein) = 2 * m * n + m

function levyarea(
        W::AbstractVector{T}, n::Integer, alg::Milstein;
        rng::AbstractRNG = default_rng()
    ) where {T <: AbstractFloat}
    m = length(W)
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    Y .= (Y .- √(T(2)) .* W) ./ (1:n)'
    A = Y * X
    # Match LevyArea.jl: M = randn!(rng, view(Y, :, 1)) — uses view into modified Y
    M = randn!(rng, view(Y, :, 1))
    a = T(√(2 * trigamma(n + 1)))
    A .+= a .* W .* M'
    for i in 1:m
        @inbounds A[i, i] = zero(T)
        for j in (i + 1):m
            @inbounds A[i, j] = (A[i, j] - A[j, i]) / (2 * T(π))
            @inbounds A[j, i] = -A[i, j]
        end
    end
    return A
end

function levyarea(
        W::AbstractVector{T}, n::Integer, alg::Milstein,
        coeffs::LevyAreaCoefficients{T}
    ) where {T <: AbstractFloat}
    M = coeffs.tail[1:(coeffs.m)]
    return _milstein_area(W, coeffs.X, copy(coeffs.Y), M, n, coeffs.m)
end

function _milstein_area(W::AbstractVector{T}, X, Y, M, n, m) where {T}
    Y = copy(Y)
    Y .= (Y .- √(T(2)) .* W) ./ (1:n)'
    A = Y * X
    a = T(√(2 * trigamma(n + 1)))
    A .+= a .* W .* M'
    for i in 1:m
        @inbounds A[i, i] = zero(T)
        for j in (i + 1):m
            @inbounds A[i, j] = (A[i, j] - A[j, i]) / (2 * T(π))
            @inbounds A[j, i] = -A[i, j]
        end
    end
    return A
end
