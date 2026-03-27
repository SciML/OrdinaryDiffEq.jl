"""
    ito_correction!(I, h=1)

Subtract h/2 from every diagonal element of I (Itô correction).
"""
function ito_correction!(I, h = 1)
    m, n = size(I)
    m == n || throw(DimensionMismatch("Matrix is not square: dimensions are $(size(I))"))
    return @inbounds for i in 1:m
        I[i, i] -= h / 2
    end
end

"""
    iterated_integrals(W::AbstractVector, h, eps=h^(3/2); ...)

Simulate iterated stochastic integrals ∫₀ʰ∫₀ˢdWᵢ(t)dWⱼ(s) for all pairs i,j.
"""
function iterated_integrals(
        W::AbstractVector{T}, h::Real, eps::Real = h^(3 / 2);
        ito_correction = true,
        error_norm::AbstractErrorNorm = MaxL2(),
        alg::AbstractIteratedIntegralAlgorithm = optimal_algorithm(length(W), h, eps, error_norm),
        rng::AbstractRNG = default_rng()
    ) where {T <: AbstractFloat}
    m = length(W)
    n = terms_needed(m, h, eps, alg, error_norm)
    I = levyarea(W / √h, n, alg; rng = rng)
    if ito_correction
        ito_correction!(I)
    end
    I .= 0.5 .* W .* W' .+ h .* I
    return I
end

"""
    iterated_integrals(W::AbstractVector, h, coeffs::LevyAreaCoefficients; ...)

Compute iterated integrals using pre-generated Fourier coefficients (deterministic).
"""
function iterated_integrals(
        W::AbstractVector{T}, h::Real,
        coeffs::LevyAreaCoefficients{T};
        ito_correction = true,
        alg::AbstractIteratedIntegralAlgorithm = MronRoe()
    ) where {T <: AbstractFloat}
    n = coeffs.n
    I = levyarea(W / √h, n, alg, coeffs)
    if ito_correction
        ito_correction!(I)
    end
    I .= 0.5 .* W .* W' .+ h .* I
    return I
end

"""
    iterated_integrals(W::AbstractVector, q_12, h, eps; ...)

Q-Wiener process version.
"""
function iterated_integrals(
        W::AbstractVector{T}, q_12::AbstractVector, h::Real, eps::Real;
        ito_correction = true,
        error_norm::AbstractErrorNorm = FrobeniusL2(),
        alg::AbstractIteratedIntegralAlgorithm = optimal_algorithm(length(W), q_12, h, eps, error_norm),
        rng::AbstractRNG = default_rng()
    ) where {T <: AbstractFloat}
    m = length(W)
    n = terms_needed(m, q_12, h, eps, alg, error_norm)
    I = levyarea(W ./ q_12 ./ √h, n, alg; rng = rng)
    if ito_correction
        ito_correction!(I)
    end
    I .= 0.5 .* W .* W' .+ h .* q_12' .* I .* q_12
    return I
end

# Scalar versions (exact)
iterated_integrals(W::Real, h::Real, eps::Real = 0.0; ito_correction = true, kwargs...) =
    ito_correction ? 0.5W^2 - 0.5h : 0.5W^2

iterated_integrals(W::Real, q_12::Real, h::Real, eps::Real; ito_correction = true, kwargs...) =
    ito_correction ? 0.5W^2 - 0.5 * h * q_12^2 : 0.5W^2
