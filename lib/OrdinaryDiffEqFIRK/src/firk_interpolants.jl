FIRK_WITH_INTERPOLATIONS = Union{RadauIIA3ConstantCache, RadauIIA3Cache, RadauIIA5ConstantCache, RadauIIA5Cache,
    RadauIIA9ConstantCache, RadauIIA9Cache, AdaptiveRadauConstantCache, AdaptiveRadauCache}

function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::FIRK_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::FIRK_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{RadauIIA3ConstantCache, RadauIIA3Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cache.cont1, cache.cont2
    @unpack cache.tab.c1
    c1m1 = c1 - 1
    @..  y₀ + Θ * (cont1 + (Θ - c1m1) * cont2)
end


@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cache.cont1, cache.cont2, cache.cont3
    @unpack cache.tab.c1, cache.tab.c2
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    @..  y₀ + Θ * (cont1 + (Θ - c2m1) * (cont2 + (Θ - c1m1) * cont3))
end

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{RadauIIA9ConstantCache, RadauIIA9Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cache.cont1, cache.cont2, cache.cont3, cache.cont4, cache.cont5
    @unpack cache.tab.c1, cache.tab.c2, cache.tab.c3, cache.tab.c4
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c3m1 = c3 - 1
    c4m1 = c4 - 1
    @..  y₀ + Θ * (cont1 + (Θ - c4m1) * (cont2 + (Θ - c3m1) * (cont3 + (Θ - c2m1) * (cont4 + (Θ - c1m1) * cont5))))
end

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{AdaptiveRadauConstantCache, AdaptiveRadauCache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cache.cont
    @unpack cache.tab.c
    tmp = cont[num_stages] * (Θ - c[1] + 1) + cont[num_stages - 1]
    j = num_stages - 2
    while j > 0
        tmp *= (Θ - c[num_stages - j] + 1) 
        tmp += cont[j]
        j = j - 1
    end
    tmp *= Θ
    @.. y₀ + tmp
end


