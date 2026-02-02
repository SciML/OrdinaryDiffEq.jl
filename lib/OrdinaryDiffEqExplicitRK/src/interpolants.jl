"""
    Interpolants for ExplicitRK methods

Implements the dense output dispatches required by OrdinaryDiffEqCore for
ExplicitRKCache and ExplicitRKConstantCache.
"""

const ExplicitRKCacheTypes = Union{ExplicitRKCache, ExplicitRKConstantCache}

# ============================================================================
# _ode_addsteps! implementations
# ============================================================================

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::ExplicitRKConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    (; A, c, stages) = cache

    if length(k) < stages || always_calc_begin
        copyat_or_push!(k, 1, f(uprev, p, t))

        for i in 2:stages
            utilde = uprev
            for j in 1:(i - 1)
                utilde = utilde + dt * A[j, i] * k[j]
            end
            copyat_or_push!(k, i, f(utilde, p, t + c[i] * dt))
        end
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::ExplicitRKCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    (; kk, tmp, tab) = cache
    (; A, c, stages) = tab

    if length(k) < stages || always_calc_begin
        f(kk[1], uprev, p, t)
        copyat_or_push!(k, 1, kk[1])

        for i in 2:stages
            @.. broadcast = false tmp = uprev
            for j in 1:(i - 1)
                @.. broadcast = false tmp = tmp + dt * A[j, i] * kk[j]
            end
            f(kk[i], tmp, p, t + c[i] * dt)
            copyat_or_push!(k, i, kk[i])
        end
    end
    nothing
end

# ============================================================================
# _ode_interpolant implementations
# ============================================================================

# Fallback for unsupported derivative orders
function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitRKCacheTypes,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitRKCacheTypes,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

@inline function get_B_interp(cache::ExplicitRKConstantCache)
    return cache.B_interp
end

@inline function get_B_interp(cache::ExplicitRKCache)
    return cache.tab.B_interp
end

@inline function get_bi(cache::ExplicitRKConstantCache)
    return cache.bi
end

@inline function get_bi(cache::ExplicitRKCache)
    return cache.tab.bi
end

# Generate interpolant methods for derivative orders 0-3
for order in 0:3
    @eval begin
        @muladd function _ode_interpolant(
                Θ, dt, y₀, y₁, k,
                cache::ExplicitRKCacheTypes,
                idxs, T::Type{Val{$order}}, differential_vars
            )
            B_interp = get_B_interp(cache)
            bi = get_bi(cache)
            return generic_rk_interpolant(Θ, dt, y₀, k, B_interp, bi; idxs = idxs, order = $order)
        end

        @muladd function _ode_interpolant!(
                out, Θ, dt, y₀, y₁, k,
                cache::ExplicitRKCacheTypes,
                idxs, T::Type{Val{$order}}, differential_vars
            )
            B_interp = get_B_interp(cache)
            bi = get_bi(cache)
            return generic_rk_interpolant!(out, Θ, dt, y₀, k, B_interp, bi; idxs = idxs, order = $order)
        end
    end
end
