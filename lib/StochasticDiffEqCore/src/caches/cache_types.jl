"""
    StochasticCompositeCache(caches, choice_function, current)

Cache of a [`StochasticCompositeAlgorithm`](@ref).

Holds one member cache per member algorithm in `caches`, the `choice_function` that
selects among them, and `current`, the index chosen for the step being taken.
[`unwrap_alg`](@ref) and [`get_current_alg_order`](@ref) read `current` so that
`perform_step!` and the adaptivity see the member that is actually running.
"""
mutable struct StochasticCompositeCache{T, F} <: StochasticDiffEqCache
    caches::T
    choice_function::F
    current::Int
end

"""
    alg_cache(
        alg, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype,
        jump_rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, f, t, dt, ::Type{Val{iip}}, verbose
    )

Construct the per-algorithm cache used by `perform_step!`.

This is the main extension point for SDE solver subpackages: each algorithm defines a
cache type (usually with the `@cache` macro, which also generates the `full_cache`,
`rand_cache`, and `ratenoise_cache` accessors used for `resize!`) and one
`alg_cache` method that allocates it. The `Val{iip}` argument selects the in-place or
out-of-place variant, and the `*NoUnits` type arguments give the element types the
error estimates should be computed in.

The fallback method defined here errors, pointing at the solver subpackage that needs
to be loaded for `alg`.
"""
function alg_cache(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}, args...)
    error("No cache constructor was loaded for $(typeof(alg)). Load the solver subpackage that defines this algorithm.")
end

function alg_cache(
        alg::algType,
        prob,
        u,
        ΔW,
        ΔZ,
        p,
        rate_prototype,
        noise_rate_prototype,
        jump_rate_prototype,
        ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits},
        uprev,
        f,
        t,
        dt,
        ::Type{Val{T}},
        verbose
    ) where {
        T, algType <: StochasticCompositeAlgorithm,
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits,
    }
    caches = map(
        (x) -> alg_cache(
            x, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype,
            jump_rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
            tTypeNoUnits, uprev, f, t, dt, Val{T}, verbose
        ),
        alg.algs
    )
    return StochasticCompositeCache(caches, alg.choice_function, 1)
end
