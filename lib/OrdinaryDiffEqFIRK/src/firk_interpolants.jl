FIRK_WITH_INTERPOLATIONS = Union{
    RadauIIA3ConstantCache, RadauIIA3Cache, RadauIIA5ConstantCache, RadauIIA5Cache,
    RadauIIA9ConstantCache, RadauIIA9Cache, AdaptiveRadauConstantCache, AdaptiveRadauCache,
}

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RadauIIA3ConstantCache, RadauIIA3Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    (; c1) = cache.tab
    z1 = k[3]
    z2 = k[4] 
    #(0,0), (c1, z1), (1,z2)
    l1 = z2 * (Θ - 0) * (Θ - c1)/ ((1 - 0) * (1 - c1))
    l2 = z1 * (Θ - 0) * (Θ - 1)/ ((c1 - 0) *  (c1 - 1))
    @.. y₀ + l1 + l2
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA3ConstantCache, RadauIIA3Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1) = cache.tab
    z1 = k[3]
    z2 = k[4] 
    #(0,0), (c1, z1), (1,z2)
    l1 = z2 * (Θ - 0) * (Θ - c1)/ ((1 - 0) * (1 - c1))
    l2 = z1 * (Θ - 0) * (Θ - 1)/ ((c1 - 0) *  (c1 - 1))
    @.. out = y₀ + l1 + l2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2) = cache.tab
    z1 = k[3]
    z2 = k[4]
    z3 = k[5]
    
    l1 = z3 * (Θ - 0) * (Θ - c1) * (Θ - c2) / ((1 - 0) * (1 - c1) * (1 - c2))
    l2 = z2 * (Θ - 0) * (Θ - c1) * (Θ - 1) / ((c2 - 0) * (c2 - c1) * (c2 - 1))
    l3 = z1 * (Θ - 0) * (Θ - c2) * (Θ - 1) / ((c1 - 0) * (c1 - c2) * (c1 - 1))
    #l4 = 0 * (Θ - c1) * (Θ - c2) * (Θ - 1) / ((0 - c1) * (0 - c2) * (0 - 1))
    @.. y₀ + l1 + l2 + l3
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA5ConstantCache, RadauIIA5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2) = cache.tab
    z1 = k[3]
    z2 = k[4]
    z3 = k[5]

    l1 = z3 * (Θ - 0) * (Θ - c1) * (Θ - c2) / ((1 - 0) * (1 - c1) * (1 - c2))
    l2 = z2 * (Θ - 0) * (Θ - c1) * (Θ - 1) / ((c2 - 0) * (c2 - c1) * (c2 - 1))
    l3 = z1 * (Θ - 0) * (Θ - c2) * (Θ - 1) / ((c1 - 0) * (c1 - c2) * (c1 - 1))
    #l4 = 0 * (Θ - c1) * (Θ - c2) * (Θ - 1) / ((0 - c1) * (0 - c2) * (0 - 1))
    @.. out = y₀ + l1 + l2 + l3
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RadauIIA9ConstantCache, RadauIIA9Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2, c3, c4) = cache.tab

    z1 = k[3]
    z2 = k[4]
    z3 = k[5]
    z4 = k[6]
    z5 = k[7]
    
    l1 = z5 * (Θ - 0) * (Θ - c1) * (Θ - c2) * (Θ - c3) * (Θ - c4) / ((1 - 0) * (1 - c1) * (1 - c2) * (1 - c3) * (1 - c4))
    l2 = z4 * (Θ - 0) * (Θ - c1) * (Θ - c2) * (Θ - c3) * (Θ - 1) / ((c4 - 0) * (c4 - c1) * (c4 - c2) * (c4 - c3) * (c4 - 1))
    l3 = z3 * (Θ - 0) * (Θ - c1) * (Θ - c2) * (Θ - c4) * (Θ - 1) / ((c3 - 0) * (c3 - c1) * (c3 - c2) * (c3 - c4) * (c3 - 1))
    l4 = z2 * (Θ - 0) * (Θ - c1) * (Θ - c3) * (Θ - c4) * (Θ - 1) / ((c2 - 0) * (c2 - c1) * (c2 - c3) * (c2 - c4) * (c2 - 1))
    l5 = z1 * (Θ - 0) * (Θ - c2) * (Θ - c3) * (Θ - c4) * (Θ - 1) / ((c1 - 0) * (c1 - c2) * (c1 - c3) * (c1 - c4) * (c1 - 1))
    #l6 = 0 * (Θ - c1) * (Θ - c2) * (Θ - c3) * (Θ - c4) * (Θ - 1) / ((0 - c1) * (0 - c2) * (0 - c3) * (0 - c4) * (0 - 1))
    @.. y₀ + l1 + l2 + l3 + l4 + l5 

end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{RadauIIA9ConstantCache, RadauIIA9Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; c1, c2, c3, c4) = cache.tab

    z1 = k[3]
    z2 = k[4]
    z3 = k[5]
    z4 = k[6]
    z5 = k[7]
    
    l1 = z5 * (Θ - 0) * (Θ - c1) * (Θ - c2) * (Θ - c3) * (Θ - c4) / ((1 - 0) * (1 - c1) * (1 - c2) * (1 - c3) * (1 - c4))
    l2 = z4 * (Θ - 0) * (Θ - c1) * (Θ - c2) * (Θ - c3) * (Θ - 1) / ((c4 - 0) * (c4 - c1) * (c4 - c2) * (c4 - c3) * (c4 - 1))
    l3 = z3 * (Θ - 0) * (Θ - c1) * (Θ - c2) * (Θ - c4) * (Θ - 1) / ((c3 - 0) * (c3 - c1) * (c3 - c2) * (c3 - c4) * (c3 - 1))
    l4 = z2 * (Θ - 0) * (Θ - c1) * (Θ - c3) * (Θ - c4) * (Θ - 1) / ((c2 - 0) * (c2 - c1) * (c2 - c3) * (c2 - c4) * (c2 - 1))
    l5 = z1 * (Θ - 0) * (Θ - c2) * (Θ - c3) * (Θ - c4) * (Θ - 1) / ((c1 - 0) * (c1 - c2) * (c1 - c3) * (c1 - c4) * (c1 - 1))

    @.. out = y₀ + l1 + l2 + l3 + l4 + l5

end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{AdaptiveRadauConstantCache, AdaptiveRadauCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; num_stages, index) = cache
    (; c) = cache.tabs[index]

    l = Vector{typeof(y₀)}(undef, num_stages)
    for i in 1:num_stages
        l[i] = k[i + 2]
        for j in 1:num_stages
            if j != i
                l[i] *= (Θ - c[j]) / (c[i] - c[j])
            end
        end
        l[i] *= Θ / c[i]
    end
    increment = sum(l)
    @.. y₀ + increment
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::Union{AdaptiveRadauConstantCache, AdaptiveRadauCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    (; num_stages, index) = cache
    (; c) = cache.tabs[index]

    l = Vector{typeof(y₀)}(undef, num_stages)
    for i in 1:num_stages
        l[i] = k[i + 2]
        for j in 1:num_stages
            if j != i
                l[i] *= (Θ - c[j]) / (c[i] - c[j])
            end
        end
        l[i] *= Θ / c[i]
    end
    increment = sum(l)
    @.. out = y₀ + increment
end
