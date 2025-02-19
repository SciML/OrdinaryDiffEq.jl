FIRK_WITH_INTERPOLATIONS = Union{RadauIIA3ConstantCache, RadauIIA3Cache, RadauIIA5ConstantCache, RadauIIA5Cache,
    RadauIIA9ConstantCache, RadauIIA9Cache, AdaptiveRadauConstantCache, AdaptiveRadauCache}

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{RadauIIA3ConstantCache, RadauIIA3Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cont1, cont2 = cache
    @unpack c1 = cache.tab
    c1m1 = c1 - 1
    @..  y₀ + Θ * (cont1 + (Θ - c1m1) * cont2)
end

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cont1, cont2, cont3 = cache
    @unpack c1, c2 = cache.tab
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    Θdt = (1-Θ)
    @.. y₁ - Θdt * (cont1 - (Θdt - c2m1) * (cont2 - (Θdt - c1m1) * cont3))
end

@muladd function _ode_interpolant!(
    out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cont1, cont2, cont3 = cache
    @unpack c1, c2 = cache.tab
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    Θdt = (1-Θ)
    @.. out = y₁ - Θdt * (cont1 - (Θdt - c2m1) * (cont2 - (Θdt - c1m1) * cont3))
end

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{RadauIIA9ConstantCache, RadauIIA9Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cont1, cont2, cont3, cont4, cont5 = cache
    @unpack c1, c2, c3, c4 = cache.tab
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c3m1 = c3 - 1
    c4m1 = c4 - 1
    @..  y₀ + Θ * (cont1 + (Θ - c4m1) * (cont2 + (Θ - c3m1) * (cont3 + (Θ - c2m1) * (cont4 + (Θ - c1m1) * cont5))))
end

@muladd function _ode_interpolant(
    Θ, dt, y₀, y₁, k, cache::Union{AdaptiveRadauConstantCache, AdaptiveRadauCache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @unpack cont = cache
    @unpack c = cache.tab
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


