"""
    generate_coefficients(m, n, alg, rng=default_rng(); T=Float64)

Generate Fourier coefficients for Lévy area computation. The layout of X and Y
depends on the algorithm:
- Fourier, Milstein, Wiktorsson: X is n×m, Y is m×n
- MronRoe: X is m×n, Y is n×m
"""
function generate_coefficients(
        m::Int, n::Int, alg::Fourier,
        rng::AbstractRNG = default_rng(); T::Type = Float64
    )
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    return LevyAreaCoefficients{T}(X, Y, T[], m, n)
end

function generate_coefficients(
        m::Int, n::Int, alg::Milstein,
        rng::AbstractRNG = default_rng(); T::Type = Float64
    )
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    M = randn(rng, T, m)
    return LevyAreaCoefficients{T}(X, Y, M, m, n)
end

function generate_coefficients(
        m::Int, n::Int, alg::Wiktorsson,
        rng::AbstractRNG = default_rng(); T::Type = Float64
    )
    X = randn(rng, T, n, m)
    Y = randn(rng, T, m, n)
    n_tail = (m^2 - m) ÷ 2
    tail = randn(rng, T, n_tail)
    return LevyAreaCoefficients{T}(X, Y, tail, m, n)
end

function generate_coefficients(
        m::Int, n::Int, alg::MronRoe,
        rng::AbstractRNG = default_rng(); T::Type = Float64
    )
    # MronRoe uses transposed layout
    X = randn(rng, T, m, n)
    Y = randn(rng, T, n, m)
    Ψ = randn(rng, T, m)
    n_tail_extra = (m^2 - m) ÷ 2
    tail_extra = randn(rng, T, n_tail_extra)
    tail = vcat(Ψ, tail_extra)
    return LevyAreaCoefficients{T}(X, Y, tail, m, n)
end
