mutable struct StochasticCompositeCache{T, F} <: StochasticDiffEqCache
    caches::T
    choice_function::F
    current::Int
end

# Concrete `alg_cache` methods live in the solver subpackages. This fallback makes
# the `alg_cache` call in `_sde_init` statically resolvable for any algorithm
# (JET) and turns a missing implementation into a descriptive error instead of a
# `MethodError`.
function alg_cache(
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm},
        prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype,
        jump_rate_prototype, ::Type, ::Type, ::Type, uprev, f, t, dt, ::Type, verbose
    )
    throw(
        ArgumentError(
            "`alg_cache` is not implemented for algorithm $(typeof(alg)). " *
                "Load the solver subpackage that implements this algorithm."
        )
    )
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
