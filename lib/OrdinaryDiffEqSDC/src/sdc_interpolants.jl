const SDC_CACHES = Union{SDCCache, SDCConstantCache}

"""
    sdc_dense_weight(dense, j, Θ, dt, order)

Weight of the node rate `zⱼ` in the `order`-th time derivative of the collocation
polynomial at `t + Θ dt`.
"""
function sdc_dense_weight(dense, j, Θ, dt, order)
    lowest = max(1, order)
    w = zero(Θ * dense[1, j])
    for i in size(dense, 1):-1:lowest
        w = w * Θ + dense[i, j] * prod((i - order + 1):i; init = 1)
    end
    return w * Θ^(lowest - order) / dt^order
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs::Nothing,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    (; dense) = cache.tab
    out = D == 0 ? y₀ : zero(y₀)
    for j in eachindex(k)
        out = out .+ sdc_dense_weight(dense, j, Θ, dt, D) .* k[j]
    end
    return out
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    (; dense) = cache.tab
    out = D == 0 ? y₀[idxs] : zero(y₀[idxs])
    for j in eachindex(k)
        out = out .+ sdc_dense_weight(dense, j, Θ, dt, D) .* k[j][idxs]
    end
    return out
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs::Nothing,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    (; dense) = cache.tab
    D == 0 ? copyto!(out, y₀) : fill!(out, false)
    for j in eachindex(k)
        w = sdc_dense_weight(dense, j, Θ, dt, D)
        kj = k[j]
        @.. broadcast = false out = out + w * kj
    end
    return out
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    (; dense) = cache.tab
    D == 0 ? copyto!(out, view(y₀, idxs)) : fill!(out, false)
    for j in eachindex(k)
        w = sdc_dense_weight(dense, j, Θ, dt, D)
        kj = view(k[j], idxs)
        @.. broadcast = false out = out + w * kj
    end
    return out
end

# `k` already holds every node rate of the step, so there is nothing to add; the
# generic method would overwrite `k[1]` and `k[2]` with endpoint derivatives.
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::SDC_CACHES, always_calc_begin = false,
        allow_calc_end = true, force_calc_end = false
    )
    return nothing
end

function SciMLBase.interp_summary(::Type{cacheType}, dense::Bool) where {cacheType <: SDC_CACHES}
    return dense ? "collocation polynomial" : "1st order linear"
end
