function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            RosenbrockCombinedConstantCache,
            RosenbrockCache,
            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
    }
    return dense ? "specialized 3rd order \"free\" stiffness-aware interpolation" :
        "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            RosenbrockCombinedConstantCache,
            RosenbrockCache,
        },
    }
    return dense ?
        "specialized 4th (Rodas6P = 5th) order \"free\" stiffness-aware interpolation" :
        "1st order linear"
end

# strip_cache for RosenbrockCache: the generic OrdinaryDiffEqCore version passes
# all-Nothing args, but several fields (dense, dtC, dtd, ks, interp_order, jac_reuse) have
# concrete types that don't accept Nothing.  Provide a custom override that clears
# only the interpolation-related fields to zero-length arrays / nothing.
function OrdinaryDiffEqCore.strip_cache(cache::RosenbrockCache)
    return RosenbrockCache(
        cache.u, cache.uprev,
        similar(cache.dense, 0),     # dense::Vector{rateType} — must not be Nothing
        cache.du, cache.du1, cache.du2,
        similar(cache.dtC, 0, 0),    # dtC::Matrix{tabType}   — must not be Nothing
        similar(cache.dtd, 0),       # dtd::Vector{tabType}   — must not be Nothing
        similar(cache.ks, 0),        # ks::Vector{rateType}   — must not be Nothing
        cache.fsalfirst, cache.fsallast, cache.dT,
        cache.J, cache.W,
        cache.tmp, cache.atmp, cache.weight,
        cache.tab,
        nothing,                     # tf
        nothing,                     # uf
        cache.linsolve_tmp, cache.linsolve,
        nothing,                     # jac_config
        nothing,                     # grad_config
        cache.reltol, cache.alg,
        cache.step_limiter!, cache.stage_limiter!,
        cache.interp_order,          # interp_order::Int — must not be Nothing
        cache.jac_reuse              # jac_reuse — preserve state across strip
    )
end
