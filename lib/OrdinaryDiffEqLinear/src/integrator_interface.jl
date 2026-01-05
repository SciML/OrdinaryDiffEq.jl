@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::LinearExponential,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp,)
end
