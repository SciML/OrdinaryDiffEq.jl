"""
    Fourier()

Basic Fourier-series approximation of Lévy areas. Its truncation error has
convergence order ``1/2`` in the number of retained Fourier terms.
"""
struct Fourier <: AbstractIteratedIntegralAlgorithm end

convorder(::Fourier) = 1 // 2
errcoeff(m, h, ::Fourier, ::MaxLp{2}) = h * √3 / (√2 * π)
norv(m, n, ::Fourier) = 2 * m * n

"""
    levyarea(W, n, alg; rng=default_rng())
    levyarea(W, n, alg, coeffs)

Approximate the antisymmetric Lévy-area matrix for the normalized Wiener
increment `W` using `n` Fourier terms and algorithm `alg`.

The keyword form draws the required random coefficients from `rng`. Passing
pre-generated [`LevyAreaCoefficients`](@ref) makes the computation
deterministic. Use [`iterated_integrals`](@ref) to obtain the complete matrix
of iterated integrals over a time step of length `h`.
"""
function levyarea(
        W::AbstractVector{T}, n::Integer, alg::Fourier;
        rng::AbstractRNG = default_rng()
    ) where {T <: AbstractFloat}
    m = length(W)
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    return _fourier_area(W, X, Y, n, m)
end

function levyarea(
        W::AbstractVector{T}, n::Integer, alg::Fourier,
        coeffs::LevyAreaCoefficients{T}
    ) where {T <: AbstractFloat}
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
