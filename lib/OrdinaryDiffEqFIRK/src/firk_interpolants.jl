FIRK_WITH_INTERPOLATIONS = Union{
    RadauIIA3ConstantCache, RadauIIA3Cache, RadauIIA5ConstantCache, RadauIIA5Cache,
    RadauIIA9ConstantCache, RadauIIA9Cache, AdaptiveRadauConstantCache, AdaptiveRadauCache,
}

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RadauIIA3ConstantCache, RadauIIA3Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; cont1, cont2) = cache
    (; c1) = cache.tab
    c1m1 = c1 - 1
    Θdt = 1 - Θ
    @.. y₁ - Θdt * (k[3] - (Θdt + c1m1) * k[4])
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA3ConstantCache, RadauIIA3Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1) = cache.tab
    c1m1 = c1 - 1
    Θdt = 1 - Θ
    @.. out = y₁ - Θdt * (k[3] - (Θdt + c1m1) * k[4])
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2) = cache.tab
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    Θdt = 1 - Θ
    @.. y₁ - Θdt * (k[3] - (Θdt + c2m1) * (k[4] - (Θdt + c1m1) * k[5]))
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2) = cache.tab
    (; dtprev) = cache
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    Θdt = 1 - Θ
    @.. out = y₁ - Θdt * (k[3] - (Θdt + c2m1) * (k[4] - (Θdt + c1m1) * k[5]))
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RadauIIA9ConstantCache, RadauIIA9Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2, c3, c4) = cache.tab
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c3m1 = c3 - 1
    c4m1 = c4 - 1
    Θdt = 1 - Θ
    @.. y₁ -
        Θdt * (
        k[3] -
            (Θdt + c4m1) *
            (k[4] - (Θdt + c3m1) * (k[5] - (Θdt + c2m1) * (k[6] - (Θdt + c1m1) * k[7])))
    )
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA9ConstantCache, RadauIIA9Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2, c3, c4) = cache.tab
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c3m1 = c3 - 1
    c4m1 = c4 - 1
    Θdt = 1 - Θ
    @.. out = y₁ -
        Θdt * (
        k[3] -
            (Θdt + c4m1) *
            (k[4] - (Θdt + c3m1) * (k[5] - (Θdt + c2m1) * (k[6] - (Θdt + c1m1) * k[7])))
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{AdaptiveRadauConstantCache, AdaptiveRadauCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; num_stages, index) = cache
    (; c) = cache.tabs[index]
    Θdt = 1 - Θ
    tmp = k[num_stages + 1] - k[num_stages + 2] * (Θdt + c[1] - 1)
    j = num_stages - 2
    while j > 0
        tmp *= (Θdt + c[num_stages - j] - 1)
        tmp = k[j + 2] - tmp
        j = j - 1
    end
    tmp *= Θdt
    @.. y₁ - tmp
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{AdaptiveRadauConstantCache, AdaptiveRadauCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; num_stages, index) = cache
    (; c) = cache.tabs[index]
    Θdt = 1 - Θ
    tmp = k[num_stages + 1] - k[num_stages + 2] * (Θdt + c[1] - 1)
    j = num_stages - 2
    while j > 0
        tmp *= (Θdt + c[num_stages - j] - 1)
        tmp = k[j + 2] - tmp
        j = j - 1
    end
    tmp *= Θdt
    @.. out = y₁ - tmp
end
