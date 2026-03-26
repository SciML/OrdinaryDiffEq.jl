struct Wiktorsson <: AbstractIteratedIntegralAlgorithm end

convorder(::Wiktorsson) = 1 // 1
errcoeff(m, h, ::Wiktorsson, ::MaxLp{2}) = √(5 * m) * h / (√12 * π)
norv(m, n, ::Wiktorsson) = 2 * m * n + (m^2 - m) ÷ 2

function levyarea(
        W::AbstractVector{T}, n::Integer, alg::Wiktorsson;
        rng::AbstractRNG = default_rng()
    ) where {T <: AbstractFloat}
    m = length(W)
    A = similar(W, m, m)
    G = similar(W, m, m)
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    Y .= (Y .- √(T(2)) .* W) ./ (1:n)'
    mul!(A, Y, X)
    a = T(√(2 * trigamma(n + 1)))
    # Match LevyArea.jl: individual randn calls in loop
    for j in 1:m
        @inbounds G[j, j] = zero(T)
        for i in (j + 1):m
            g = a * randn(rng, T)
            @inbounds G[i, j] = g
            @inbounds G[j, i] = -g
            @inbounds A[i, j] += g
        end
    end
    A .+= inv(1 + √(1 + W' * W)) .* (G * W) .* W'
    G .= inv(2 * T(π)) .* (A .- A')
    return G
end

function levyarea(
        W::AbstractVector{T}, n::Integer, alg::Wiktorsson,
        coeffs::LevyAreaCoefficients{T}
    ) where {T <: AbstractFloat}
    return _wiktorsson_area(W, coeffs.X, copy(coeffs.Y), coeffs.tail, n, coeffs.m)
end

function _wiktorsson_area(W::AbstractVector{T}, X, Y, tail, n, m) where {T}
    Y = copy(Y)
    A = similar(W, m, m)
    G = similar(W, m, m)
    Y .= (Y .- √(T(2)) .* W) ./ (1:n)'
    mul!(A, Y, X)
    a = T(√(2 * trigamma(n + 1)))
    idx = 0
    for j in 1:m
        @inbounds G[j, j] = zero(T)
        for i in (j + 1):m
            idx += 1
            g = a * tail[idx]
            @inbounds G[i, j] = g
            @inbounds G[j, i] = -g
            @inbounds A[i, j] += g
        end
    end
    A .+= inv(1 + √(1 + W' * W)) .* (G * W) .* W'
    G .= inv(2 * T(π)) .* (A .- A')
    return G
end
