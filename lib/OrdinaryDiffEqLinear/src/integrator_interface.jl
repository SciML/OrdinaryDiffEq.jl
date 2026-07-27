# Core's `OrdinaryDiffEqLinearExponentialAlgorithm` overload returns
# `(cache.tmp,)`, but the migrated Linear caches expose their state scratch
# through `tmp_cache`. Specialize on the migrated cache supertype
# (`LinearMutableCache` — strictly more specific than Core's
# `OrdinaryDiffEqMutableCache`, so no ambiguity) to cover the whole
# Magnus/Lie/CG/RKMK/LinearExponential family in one method.
@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqLinearExponentialAlgorithm,
        cache::LinearMutableCache
    )
    return (cache.tmp_cache.tmp,)
end
