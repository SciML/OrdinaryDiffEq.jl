# One `NLSolver` per implicit node rather than one shared solver, because each
# node has its own `γ_m = QΔ[m,m]` and hence its own `W`. Sharing would be wrong,
# not merely slow: `do_newJW` decides `W` reuse from `J_t` without comparing `γ`
# against the one the stored `W` was built with, so every node after the first at
# a given `t` would silently reuse the first node's `W`. `OrdinaryDiffEqPDIRK`
# builds a solver vector for the same reason.

struct SDCConstantCache{N, TabType} <: OrdinaryDiffEqConstantCache
    nlsolvers::N
    tab::TabType
    solver_index::Vector{Int}
end

@cache mutable struct SDCCache{uType, rateType, uNoUnitsType, N, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    ulow::uType
    atmp::uNoUnitsType
    # One entry per node so the sweep can run the nodes concurrently.
    tmp::Vector{uType}
    ubuf::Vector{uType}
    k::Vector{rateType}
    z::Vector{uType}
    znew::Vector{uType}
    # Per-node solve outcome, collected instead of returning from inside the loop.
    failed::Vector{Bool}
    nf::Vector{Int}
    nlsolvers::N
    tab::TabType
    solver_index::Vector{Int}
end

get_fsalfirstlast(cache::SDCCache, u) = (nothing, nothing)

"""
    sdc_solver_index(QΔ)

Map each node to its position in the solver vector, or to `0` when the node is
explicit because `QΔ[m, m]` vanishes. That happens for every node of a strictly
lower triangular sweeper, and for the first node whenever `τ₁ = 0`.
"""
function sdc_solver_index(QΔ::AbstractVector{<:AbstractMatrix})
    M = size(first(QΔ), 1)
    index = zeros(Int, M)
    implicit = 0
    for m in 1:M
        # A sweep-dependent preconditioner only needs a solver for node m if some
        # sweep actually makes it implicit.
        if any(Q -> !iszero(Q[m, m]), QΔ)
            implicit += 1
            index[m] = implicit
        end
    end
    return index
end

function alg_cache(
        alg::SDC, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = SDCTableau(
        constvalue(uBottomEltypeNoUnits), alg.num_nodes, alg.node_type,
        alg.quad_type, alg.sweeper, alg.num_sweeps
    )
    M = alg.num_nodes
    solver_index = sdc_solver_index(tab.QΔ)
    nlsolvers = [
        build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
            uBottomEltypeNoUnits, tTypeNoUnits, first(tab.QΔ)[m, m], tab.nodes[m],
            Val(true), verbose
        ) for m in 1:M if !iszero(solver_index[m])
    ]
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    return SDCCache(
        u, uprev, zero(u), atmp,
        [zero(u) for _ in 1:M], [zero(u) for _ in 1:M],
        [zero(rate_prototype) for _ in 1:M],
        [zero(u) for _ in 1:M], [zero(u) for _ in 1:M],
        fill(false, M), zeros(Int, M),
        nlsolvers, tab, solver_index
    )
end

function alg_cache(
        alg::SDC, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = SDCTableau(
        constvalue(uBottomEltypeNoUnits), alg.num_nodes, alg.node_type,
        alg.quad_type, alg.sweeper, alg.num_sweeps
    )
    solver_index = sdc_solver_index(tab.QΔ)
    nlsolvers = [
        build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
            uBottomEltypeNoUnits, tTypeNoUnits, first(tab.QΔ)[m, m], tab.nodes[m],
            Val(false), verbose
        ) for m in 1:(alg.num_nodes) if !iszero(solver_index[m])
    ]
    return SDCConstantCache(nlsolvers, tab, solver_index)
end

# Needed because the generic `OrdinaryDiffEqNewtonAdaptiveAlgorithm` fallback
# returns `(cache.nlsolver.tmp, cache.atmp)`, which assumes a single solver.
SciMLBase.get_tmp_cache(integrator, ::SDC, cache::SDCCache) =
    (first(cache.tmp), first(cache.ubuf))
