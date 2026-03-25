"""
    LevyAreaCoefficients{T}

Stores the pre-generated random numbers (Fourier coefficients) needed by a
Lévy area algorithm. Separating coefficient generation from area computation
enables:
- Deterministic Lévy area given the same coefficients
- Path reconstruction from the Fourier expansion
- Consistent iterated integrals across step-size changes
"""
struct LevyAreaCoefficients{T <: AbstractFloat}
    X::Matrix{T}       # Fourier coefficients (layout depends on algorithm)
    Y::Matrix{T}       # Fourier coefficients (layout depends on algorithm)
    tail::Vector{T}    # Tail approximation randoms (algorithm-dependent)
    m::Int             # noise dimension
    n::Int             # truncation level
end

"""
    coefficient_length(m, n, alg)

Total number of random values stored in `LevyAreaCoefficients` for the given algorithm.
Equals `norv(m, n, alg)`.
"""
coefficient_length(m, n, alg::AbstractIteratedIntegralAlgorithm) = norv(m, n, alg)
