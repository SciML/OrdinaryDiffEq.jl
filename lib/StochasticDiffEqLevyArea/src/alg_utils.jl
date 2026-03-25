"""
    convorder(alg::AbstractIteratedIntegralAlgorithm)

Returns the convergence order of the algorithm w.r.t. the truncation parameter.
"""
function convorder end

"""
    errcoeff(dim, stepsize, alg, norm)

Returns the error coefficient for the truncation parameter bound.
"""
function errcoeff end

"""
    norv(dim, n, alg::AbstractIteratedIntegralAlgorithm)

Returns the number of random numbers needed to simulate the iterated integrals
for a Wiener process of dimension `dim` and with truncation parameter `n`.
"""
function norv end

"""
    terms_needed(dim, stepsize, eps, alg, norm)

Returns the number of terms needed to achieve error at most `eps`.
"""
function terms_needed(dim, stepsize, eps, alg::AbstractIteratedIntegralAlgorithm, norm::AbstractErrorNorm)
    ceil(Int64, (errcoeff(dim, stepsize, alg, norm) / eps)^(1 // convorder(alg)))
end

function terms_needed(dim, q_12, stepsize, eps, alg::AbstractIteratedIntegralAlgorithm, norm::AbstractErrorNorm)
    length(q_12) == dim || throw(ArgumentError("Length of q_12 must be equal to the dimension."))
    ceil(Int64, (errcoeff(dim, q_12, stepsize, alg, norm) / eps)^(1 // convorder(alg)))
end

"""
    effective_cost(dim, stepsize, eps, alg, norm)

Returns the number of random numbers needed with the given parameters.
"""
function effective_cost(dim, stepsize, eps, alg, norm)
    norv(dim, terms_needed(dim, stepsize, eps, alg, norm), alg)
end
function effective_cost(dim, q_12, stepsize, eps, alg, norm)
    norv(dim, terms_needed(dim, q_12, stepsize, eps, alg, norm), alg)
end

"""
    optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm=MaxL2())

Returns the optimal algorithm (fewest random numbers) for the given parameters.
"""
function optimal_algorithm(dim, stepsize, eps = stepsize^(3 / 2), norm::AbstractErrorNorm = MaxL2())
    ind = argmin([effective_cost(dim, stepsize, eps, alg, norm) for alg in ITER_INT_ALGS])
    return ITER_INT_ALGS[ind]
end
function optimal_algorithm(dim, q_12, stepsize, eps, norm::AbstractErrorNorm = FrobeniusL2())
    ind = argmin([effective_cost(dim, q_12, stepsize, eps, alg, norm) for alg in ITER_INT_ALGS])
    return ITER_INT_ALGS[ind]
end
