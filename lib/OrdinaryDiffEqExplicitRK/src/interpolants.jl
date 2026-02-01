"""
    Interpolants for ExplicitRK methods

Implements the dense output dispatches required by OrdinaryDiffEqCore for
ExplicitRKCache and ExplicitRKConstantCache.

The interpolation uses the stage derivatives from `k` (populated by _ode_addsteps!)
and the interpolation coefficients from `B_interp` in the tableau.
"""

# Union type for both cache types
const ExplicitRKCacheTypes = Union{ExplicitRKCache, ExplicitRKConstantCache}

# ============================================================================
# _ode_addsteps! implementations
# ============================================================================
# These functions recompute all stage derivatives when dense output (interpolation)
# is requested. This is necessary because during perform_step!, only k[1] and k[2]
# (first and last derivatives) are saved to integrator.k, but RK interpolation
# requires all stage derivatives.

"""
    _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::ExplicitRKConstantCache, ...)

Out-of-place version: Recomputes all stage derivatives for the ExplicitRK method
and stores them in k for use by the interpolant.
"""
@muladd function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::ExplicitRKConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    @unpack A, c, stages = cache

    if length(k) < stages || always_calc_begin
        # Compute first stage derivative
        copyat_or_push!(k, 1, f(uprev, p, t))

        # Compute middle stages
        for i in 2:stages
            # Accumulate: utilde = uprev + dt * Σⱼ A[j,i] * k[j]
            utilde = uprev
            for j in 1:(i - 1)
                utilde = utilde + dt * A[j, i] * k[j]
            end
            copyat_or_push!(k, i, f(utilde, p, t + c[i] * dt))
        end
    end
    nothing
end

"""
    _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::ExplicitRKCache, ...)

In-place version: Recomputes all stage derivatives for the ExplicitRK method
and stores them in k for use by the interpolant.
"""
@muladd function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::ExplicitRKCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)
    @unpack kk, tmp, tab = cache
    @unpack A, c, stages = tab

    if length(k) < stages || always_calc_begin
        # Compute first stage derivative (in-place)
        f(kk[1], uprev, p, t)
        copyat_or_push!(k, 1, kk[1])

        # Compute middle stages
        for i in 2:stages
            # Accumulate: tmp = uprev + dt * Σⱼ A[j,i] * kk[j]
            @.. broadcast=false tmp = uprev
            for j in 1:(i - 1)
                @.. broadcast=false tmp = tmp + dt * A[j, i] * kk[j]
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
function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::ExplicitRKCacheTypes,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::ExplicitRKCacheTypes,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

# Helper to get B_interp from either cache type
@inline function get_B_interp(cache::ExplicitRKConstantCache)
    cache.B_interp
end

@inline function get_B_interp(cache::ExplicitRKCache)
    cache.tab.B_interp
end

# Generate interpolant methods for derivative orders 0-3
for order in 0:3
    @eval begin
        # Out-of-place interpolation
        @muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                cache::ExplicitRKCacheTypes,
                idxs, T::Type{Val{$order}}, differential_vars)
            B_interp = get_B_interp(cache)
            # Use k (from integrator.k, populated by _ode_addsteps!) for interpolation
            return generic_rk_interpolant(Θ, dt, y₀, k, B_interp; idxs=idxs, order=$order)
        end

        # In-place interpolation
        @muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                cache::ExplicitRKCacheTypes,
                idxs, T::Type{Val{$order}}, differential_vars)
            B_interp = get_B_interp(cache)
            # Use k (from integrator.k, populated by _ode_addsteps!) for interpolation
            return generic_rk_interpolant!(out, Θ, dt, y₀, k, B_interp; idxs=idxs, order=$order)
        end
    end
end
