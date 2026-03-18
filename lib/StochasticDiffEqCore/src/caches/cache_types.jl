mutable struct StochasticCompositeCache{T, F} <: StochasticDiffEqCache
    caches::T
    choice_function::F
    current::Int
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
