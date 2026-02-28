"""
    Dense output interpolants for ExplicitTaylor and ExplicitTaylor2 methods.

The Taylor polynomial computed at each step provides a natural high-order
interpolant. For a step from t_n to t_n + dt with Taylor coefficients
c₁, c₂, ..., c_P stored in k[1], ..., k[P], the interpolant at t_n + h is:

    u(t_n + h) = y₀ + c₁·h + c₂·h² + ⋯ + c_P·h^P

where h = Θ·dt and Θ ∈ [0,1].
"""

const ExplicitTaylor2Caches = Union{ExplicitTaylor2Cache, ExplicitTaylor2ConstantCache}
const ExplicitTaylorCaches = Union{ExplicitTaylorCache, ExplicitTaylorConstantCache}
const AllTaylorCaches = Union{ExplicitTaylor2Caches, ExplicitTaylorCaches}

# ============================================================================
# interp_summary
# ============================================================================

function SciMLBase.interp_summary(
        ::Type{cacheType}, dense::Bool
    ) where {cacheType <: ExplicitTaylor2Caches}
    return dense ? "2nd order Taylor polynomial" : "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType}, dense::Bool
    ) where {cacheType <: ExplicitTaylorCaches}
    return dense ? "Taylor polynomial" : "1st order linear"
end

# ============================================================================
# _ode_addsteps! — no-op since coefficients are saved during perform_step!
# ============================================================================

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::AllTaylorCaches,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    return nothing
end

# ============================================================================
# ExplicitTaylor2 interpolation (order 2)
#
# k[1] = first derivative, k[2] = second derivative
# u(h) = y₀ + k[1]·h + (k[2]/2)·h²
# ============================================================================

# --- Val{0}: function value ---

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    h = Θ * dt
    @inbounds @.. y₀ + k[1] * h + k[2] / 2 * h^2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs, T::Type{Val{0}}, differential_vars
    )
    h = Θ * dt
    @inbounds @views @.. y₀[idxs] + k[1][idxs] * h + k[2][idxs] / 2 * h^2
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    h = Θ * dt
    @inbounds @.. out = y₀ + k[1] * h + k[2] / 2 * h^2
    return out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs, T::Type{Val{0}}, differential_vars
    )
    h = Θ * dt
    @inbounds @views @.. out = y₀[idxs] + k[1][idxs] * h + k[2][idxs] / 2 * h^2
    return out
end

# --- Val{1}: first derivative ---

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    h = Θ * dt
    @inbounds @.. k[1] + k[2] * h
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs, T::Type{Val{1}}, differential_vars
    )
    h = Θ * dt
    @inbounds @views @.. k[1][idxs] + k[2][idxs] * h
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    h = Θ * dt
    @inbounds @.. out = k[1] + k[2] * h
    return out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs, T::Type{Val{1}}, differential_vars
    )
    h = Θ * dt
    @inbounds @views @.. out = k[1][idxs] + k[2][idxs] * h
    return out
end

# --- Higher derivatives: throw error ---

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylor2Caches,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

# ============================================================================
# ExplicitTaylor{P} interpolation (arbitrary order)
#
# k[1] = c₁, k[2] = c₂, ..., k[P] = c_P (Taylor polynomial coefficients)
# u(h) = y₀ + c₁·h + c₂·h² + ⋯ + c_P·h^P
#
# For interpolation, use Horner's method:
# u(h) = y₀ + h·(c₁ + h·(c₂ + h·(⋯ + h·c_P)))
# ============================================================================

# --- Val{0}: function value using Horner's method ---

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorConstantCache{P},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    ) where {P}
    h = Θ * dt
    # Horner's method: y₀ + h*(k[1] + h*(k[2] + ... + h*k[P]))
    @inbounds begin
        result = k[P]
        for i in (P - 1):-1:1
            result = @.. k[i] + h * result
        end
        return @.. y₀ + h * result
    end
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorConstantCache{P},
        idxs, T::Type{Val{0}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = k[P][idxs]
        for i in (P - 1):-1:1
            result = @.. k[i][idxs] + h * result
        end
        return @.. y₀[idxs] + h * result
    end
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCache{P},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = copy(k[P])
        for i in (P - 1):-1:1
            @.. result = k[i] + h * result
        end
        return @.. y₀ + h * result
    end
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCache{P},
        idxs, T::Type{Val{0}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = k[P][idxs]
        for i in (P - 1):-1:1
            result = @.. k[i][idxs] + h * result
        end
        return @.. y₀[idxs] + h * result
    end
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCaches,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    P = _taylor_order(cache)
    h = Θ * dt
    @inbounds begin
        @.. out = k[P]
        for i in (P - 1):-1:1
            @.. out = k[i] + h * out
        end
        @.. out = y₀ + h * out
    end
    return out
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCaches,
        idxs, T::Type{Val{0}}, differential_vars
    )
    P = _taylor_order(cache)
    h = Θ * dt
    @inbounds @views begin
        @.. out = k[P][idxs]
        for i in (P - 1):-1:1
            @.. out = k[i][idxs] + h * out
        end
        @.. out = y₀[idxs] + h * out
    end
    return out
end

# --- Val{1}: first derivative ---
# du/dt = c₁ + 2c₂·h + 3c₃·h² + ⋯ + P·c_P·h^(P-1)
# Horner form: c₁ + h·(2c₂ + h·(3c₃ + ⋯ + h·P·c_P))

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorConstantCache{P},
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = @.. P * k[P]
        for i in (P - 1):-1:1
            result = @.. i * k[i] + h * result
        end
        return result
    end
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorConstantCache{P},
        idxs, T::Type{Val{1}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = @.. P * k[P][idxs]
        for i in (P - 1):-1:1
            result = @.. i * k[i][idxs] + h * result
        end
        return result
    end
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCache{P},
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = @.. P * k[P]
        for i in (P - 1):-1:1
            @.. result = i * k[i] + h * result
        end
        return result
    end
end

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCache{P},
        idxs, T::Type{Val{1}}, differential_vars
    ) where {P}
    h = Θ * dt
    @inbounds begin
        result = @.. P * k[P][idxs]
        for i in (P - 1):-1:1
            result = @.. i * k[i][idxs] + h * result
        end
        return result
    end
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCaches,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    P = _taylor_order(cache)
    h = Θ * dt
    @inbounds begin
        @.. out = P * k[P]
        for i in (P - 1):-1:1
            @.. out = i * k[i] + h * out
        end
    end
    return out
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCaches,
        idxs, T::Type{Val{1}}, differential_vars
    )
    P = _taylor_order(cache)
    h = Θ * dt
    @inbounds @views begin
        @.. out = P * k[P][idxs]
        for i in (P - 1):-1:1
            @.. out = i * k[i][idxs] + h * out
        end
    end
    return out
end

# --- Higher derivatives: throw error for unsupported orders ---

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCaches,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ExplicitTaylorCaches,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

# ============================================================================
# Utility
# ============================================================================

_taylor_order(::ExplicitTaylorCache{P}) where {P} = P
_taylor_order(::ExplicitTaylorConstantCache{P}) where {P} = P
