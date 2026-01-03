"""
    Interpolants for ExplicitRK methods

Implements the dense output dispatches required by OrdinaryDiffEqCore for
ExplicitRKCache and ExplicitRKConstantCache.
"""

# Union type for both cache types
const ExplicitRKCacheTypes = Union{ExplicitRKCache, ExplicitRKConstantCache}

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
@inline function get_B_interp(cache)
    hasproperty(cache, :B_interp) ? cache.B_interp : cache.tab.B_interp
end

# Generate interpolant methods for derivative orders 0-3
for order in 0:3
    @eval begin
        # Out-of-place interpolation
        @muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                cache::ExplicitRKCacheTypes,
                idxs, T::Type{Val{$order}}, differential_vars)
            B_interp = get_B_interp(cache)
            return generic_rk_interpolant(Θ, dt, y₀, k, B_interp; idxs=idxs, order=$order)
        end

        # In-place interpolation
        @muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                cache::ExplicitRKCacheTypes,
                idxs, T::Type{Val{$order}}, differential_vars)
            B_interp = get_B_interp(cache)
            return generic_rk_interpolant!(out, Θ, dt, y₀, k, B_interp; idxs=idxs, order=$order)
        end
    end
end