"""
    abstract type AbstractErrorNorm end

Abstract type for different error norms for iterated integral approximation.
"""
abstract type AbstractErrorNorm end

# Maximum of entry-wise Lp norms
struct MaxLp{p} <: AbstractErrorNorm end

function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::MaxLp{p}) where {p}
    maxqq = maximum(q_12[i] * q_12[j] for i in 1:m for j in 1:(i - 1))
    return maxqq * errcoeff(m, stepsize, alg, MaxLp{p}())
end

# Schatten-q norm of matrix of entry-wise Lp norms
struct SchattenqLp{p, q} <: AbstractErrorNorm end

function errcoeff(m, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::SchattenqLp{p, q}) where {p, q}
    return (m^2 - m)^(1 / q) * errcoeff(m, stepsize, alg, MaxLp{p}())
end
function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::SchattenqLp{p, q}) where {p, q}
    trQq_sq = abs2(sum(x -> x^q, q_12))
    tr_Qqsq = sum(x -> x^(2q), q_12)
    return (trQq_sq - tr_Qqsq)^(1 / q) * errcoeff(m, stepsize, alg, MaxLp{p}())
end

# Lp norm of maximum of matrix
struct LpMax{p} <: AbstractErrorNorm end

function errcoeff(m, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::LpMax{p}) where {p}
    return 1 / (2^(1 / p)) * errcoeff(m, stepsize, alg, SchattenqLp{p, p}())
end
function errcoeff(m, q_12, stepsize, alg::AbstractIteratedIntegralAlgorithm, ::LpMax{p}) where {p}
    return 1 / (2^(1 / p)) * errcoeff(m, q_12, stepsize, alg, SchattenqLp{p, p}())
end

const MaxL2 = MaxLp{2}
const FrobeniusL2 = SchattenqLp{2, 2}
