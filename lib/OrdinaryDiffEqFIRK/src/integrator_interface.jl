@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::Union{RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau, GaussLegendre},
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.atmp)
end
