@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::LinearExponential,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp,)
end