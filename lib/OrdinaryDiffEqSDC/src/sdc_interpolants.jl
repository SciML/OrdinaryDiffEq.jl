const SDC_CACHES = Union{SDCCache, SDCConstantCache}

# `k` carries the node rates only while they describe the step being interpolated; anything
# that invalidates them leaves the `k[1]`, `k[2]` pair behind instead.
sdc_has_nodes(cache, k) = length(k) == length(cache.tab.nodes) + 2

"""
    sdc_basis_value(dense, j, Θ)

The Lagrange basis polynomial `ℓⱼ(Θ)`, from the integral coefficients in `dense`.
"""
function sdc_basis_value(dense, j, Θ)
    w = zero(Θ * dense[1, j])
    for i in size(dense, 1):-1:1
        w = w * Θ + dense[i, j] * i
    end
    return w
end

"""
    sdc_dense_weight(dense, j, Θ, dt, order)

Weight of the node rate `fⱼ` in the `order`-th time derivative of the collocation
polynomial at `t + Θ dt`.
"""
function sdc_dense_weight(dense, j, Θ, dt, order)
    lowest = max(1, order)
    w = zero(Θ * dense[1, j])
    for i in size(dense, 1):-1:lowest
        w = w * Θ + dense[i, j] * prod((i - order + 1):i; init = 1)
    end
    return w * Θ^(lowest - order) * dt^(1 - order)
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs::Nothing,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    sdc_has_nodes(cache, k) || return hermite_interpolant(
        Θ, dt, y₀, y₁, k, Val{false}, idxs, T,
        interpolation_differential_vars(differential_vars, y₀, idxs)
    )
    (; dense) = cache.tab
    out = D == 0 ? y₀ : zero(y₀)
    for j in 1:size(dense, 2)
        out = out .+ sdc_dense_weight(dense, j, Θ, dt, D) .* k[j + 2]
    end
    return out
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    sdc_has_nodes(cache, k) || return hermite_interpolant(
        Θ, dt, y₀, y₁, k, Val{false}, idxs, T,
        interpolation_differential_vars(differential_vars, y₀, idxs)
    )
    (; dense) = cache.tab
    out = D == 0 ? y₀[idxs] : zero(y₀[idxs])
    for j in 1:size(dense, 2)
        out = out .+ sdc_dense_weight(dense, j, Θ, dt, D) .* k[j + 2][idxs]
    end
    return out
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs::Nothing,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    sdc_has_nodes(cache, k) || return hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs, T,
        interpolation_differential_vars(differential_vars, y₀, idxs)
    )
    (; dense) = cache.tab
    D == 0 ? copyto!(out, y₀) : fill!(out, false)
    for j in 1:size(dense, 2)
        w = sdc_dense_weight(dense, j, Θ, dt, D)
        kj = k[j + 2]
        @.. broadcast = false out = out + w * kj
    end
    return out
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::SDC_CACHES, idxs,
        T::Type{Val{D}}, differential_vars
    ) where {D}
    sdc_has_nodes(cache, k) || return hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs, T,
        interpolation_differential_vars(differential_vars, y₀, idxs)
    )
    (; dense) = cache.tab
    D == 0 ? copyto!(out, view(y₀, idxs)) : fill!(out, false)
    for j in 1:size(dense, 2)
        w = sdc_dense_weight(dense, j, Θ, dt, D)
        kj = view(k[j + 2], idxs)
        @.. broadcast = false out = out + w * kj
    end
    return out
end

# A recompute request means the step changed under the node rates (a callback shortened it),
# so they are dropped for the endpoint pair the generic interpolant uses.
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::SDC_CACHES, always_calc_begin = false,
        allow_calc_end = true, force_calc_end = false
    )
    (always_calc_begin || length(k) < 2) || return nothing
    if cache isa OrdinaryDiffEqMutableCache
        rtmp = similar(u, eltype(eltype(k)))
        f(rtmp, uprev, p, t)
        copyat_or_push!(k, 1, rtmp)
        f(rtmp, u, p, t + dt)
        copyat_or_push!(k, 2, rtmp)
    else
        copyat_or_push!(k, 1, f(uprev, p, t))
        copyat_or_push!(k, 2, f(u, p, t + dt))
    end
    resize!(k, 2)
    return nothing
end

function SciMLBase.interp_summary(::Type{cacheType}, dense::Bool) where {cacheType <: SDC_CACHES}
    return dense ? "collocation polynomial" : "1st order linear"
end
