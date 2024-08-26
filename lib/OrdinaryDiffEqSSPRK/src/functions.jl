@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::Union{SSPRK22, SSPRK33, SSPRK53_2N1,
            SSPRK53_2N2, SSPRK43, SSPRK432,
            SSPRK932},
        cache::OrdinaryDiffEqMutableCache)
    (cache.k,)
end
