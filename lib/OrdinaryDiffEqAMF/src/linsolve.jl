struct SciMLOpFactorization <: LinearSolve.SciMLLinearSolveAlgorithm end

function LinearSolve.init_cacheval(
    alg::SciMLOpFactorization,
    A::SciMLOperators.AbstractSciMLOperator,
    b,
    u,
    Pl,
    Pr,
    maxiters::Int,
    abstol,
    reltol,
    verbose::Union{LinearVerbosity, Bool},
    assumptions::LinearSolve.OperatorAssumptions,
)
    _fact = LinearAlgebra.factorize(A)
    return cache_operator(_fact, u)
end

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::SciMLOpFactorization; kwargs...)
    if cache.isfresh
        _fact = LinearAlgebra.factorize(cache.A)
        cache.cacheval = cache_operator(_fact, cache.u)
        cache.isfresh = false
    end
    y = ldiv!(cache.u, cache.cacheval, cache.b)
    return SciMLBase.build_linear_solution(alg, y, nothing, cache)
end

LinearSolve.needs_concrete_A(::SciMLOpFactorization) = true

default_amf_linsolve() = SciMLOpFactorization()
