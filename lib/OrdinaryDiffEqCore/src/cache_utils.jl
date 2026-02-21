is_constant_cache(::OrdinaryDiffEqConstantCache) = true
is_constant_cache(::OrdinaryDiffEqCache) = false
is_constant_cache(cache::CompositeCache) = is_constant_cache(cache.caches[1])
function is_constant_cache(::DefaultCache{Cache1}) where {Cache1}
    return Cache1 <: OrdinaryDiffEqConstantCache
end
function is_constant_cache(::DefaultCacheVF64{Cache1}) where {Cache1}
    return Cache1 <: OrdinaryDiffEqConstantCache
end

function SciMLBase.unwrap_cache(integrator::ODEIntegratorType, is_stiff)
    alg = integrator.alg
    cache = integrator.cache
    iscomp = alg isa CompositeAlgorithm
    if !iscomp
        return cache
    elseif cache isa DefaultCacheType
        current = integrator.cache.current
        if current == 1
            return cache.cache1
        elseif current == 2
            return cache.cache2
        elseif current == 3
            return cache.cache3
        elseif current == 4
            return cache.cache4
        elseif current == 5
            return cache.cache5
        elseif current == 6
            return cache.cache6
        else
            error("This should not occur (please report a bug)")
        end
    elseif alg.choice_function isa AutoSwitch
        num = is_stiff ? 2 : 1
        return cache.caches[num]
    else
        return cache.caches[integrator.cache.current]
    end
end

@deprecate alg_cache(
    alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, u, rate_prototype,
    uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f,
    t, dt, reltol, p, calck, ::Type{Val{iip}}
) where {iip} alg_cache(
    alg,
    u,
    rate_prototype,
    uEltypeNoUnits,
    uBottomEltypeNoUnits,
    tTypeNoUnits,
    uprev,
    uprev2,
    f, t,
    dt,
    reltol,
    p,
    calck,
    Val(iip)
)
