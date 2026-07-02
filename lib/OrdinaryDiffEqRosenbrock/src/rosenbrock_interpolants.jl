### Fallbacks to capture
ROSENBROCKS_WITH_INTERPOLATIONS = Union{
    Rosenbrock23ConstantCache, Rosenbrock23Cache,
    Rosenbrock32ConstantCache, Rosenbrock32Cache,
    RosenbrockCombinedConstantCache,
    RosenbrockCache,
    HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
}

# Mask used by the `interp_order == -1` Hermite fallback (methods with no
# H matrix store k[1]=f₀, k[2]=f₁). For DAE algebraic variables, f₀/f₁ are
# residuals, not derivatives, so the Hermite correction is invalid there —
# the mask is zero on algebraic vars to fall back to linear interpolation.
@inline _dae_hermite_mask(::Nothing, ::Any) = true
@inline _dae_hermite_mask(::DifferentialVarsUndefined, ::Any) = false
@inline _dae_hermite_mask(dv::AbstractArray, ::Nothing) = dv
@inline _dae_hermite_mask(dv::AbstractArray, idxs::Number) = dv[idxs]
@inline _dae_hermite_mask(dv::AbstractArray, idxs) = @view dv[idxs]

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::ROSENBROCKS_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::ROSENBROCKS_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

####

"""
From MATLAB ODE Suite by Shampine
"""
@def rosenbrock2332unpack begin
    if cache isa OrdinaryDiffEqMutableCache
        d = cache.tab.d
    else
        d = cache.d
    end
end

@def rosenbrock2332pre0 begin
    @rosenbrock2332unpack
    c1 = Θ * (1 - Θ) / (1 - 2d)
    c2 = Θ * (Θ - 2d) / (1 - 2d)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache,
        }, idxs::Nothing,
        T::Type{Val{0}}, differential_vars
    )
    @rosenbrock2332pre0
    @inbounds y₀ + dt * (c1 * k[1] + c2 * k[2])
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23Cache, Rosenbrock32Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    @rosenbrock2332pre0
    @inbounds @.. y₀ + dt * (c1 * k[1] + c2 * k[2])
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache, Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs, T::Type{Val{0}}, differential_vars
    )
    @rosenbrock2332pre0
    @.. y₀[idxs] + dt * (c1 * k[1][idxs] + c2 * k[2][idxs])
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    @rosenbrock2332pre0
    @inbounds @.. out = y₀ + dt * (c1 * k[1] + c2 * k[2])
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs, T::Type{Val{0}}, differential_vars
    )
    @rosenbrock2332pre0
    @views @.. out = y₀[idxs] + dt * (c1 * k[1][idxs] + c2 * k[2][idxs])
    out
end

# First Derivative of the dense output
@def rosenbrock2332pre1 begin
    @rosenbrock2332unpack
    c1diff = (1 - 2 * Θ) / (1 - 2 * d)
    c2diff = (2 * Θ - 2 * d) / (1 - 2 * d)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache, Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    @rosenbrock2332pre1
    @.. c1diff * k[1] + c2diff * k[2]
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache, Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs, T::Type{Val{1}}, differential_vars
    )
    @rosenbrock2332pre1
    @.. c1diff * k[1][idxs] + c2diff * k[2][idxs]
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    @rosenbrock2332pre1
    @.. out = c1diff * k[1] + c2diff * k[2]
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache,
        }, idxs, T::Type{Val{1}}, differential_vars
    )
    @rosenbrock2332pre1
    @views @.. out = c1diff * k[1][idxs] + c2diff * k[2][idxs]
    out
end

"""
From MATLAB ODE Suite by Shampine
"""

@muladd function _ode_interpolant(
        Θ, dt,
        y₀,
        y₁,
        k,
        cache::Union{
            RosenbrockCombinedConstantCache,
            RosenbrockCache,
            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    Θ1 = 1 - Θ
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        # Generic Hermite: k[1]=f₀, k[2]=f₁ (methods with no H matrix).
        # Mask the Hermite correction by `differential_vars` so DAE algebraic
        # variables fall back to linear interpolation (f₀/f₁ are residuals there).
        m = _dae_hermite_mask(differential_vars, nothing)
        if m === true
            @.. Θ1 * y₀ + Θ * y₁ +
                Θ * (Θ - 1) * ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
        elseif m === false
            @.. Θ1 * y₀ + Θ * y₁
        else
            @.. Θ1 * y₀ + Θ * y₁ +
                m * Θ * (Θ - 1) * ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
    elseif cache.interp_order == 4
        @.. Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * (k[3] + Θ * k[4]))))
    else
        @.. Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
    end
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            RosenbrockCombinedConstantCache, RosenbrockCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs, T::Type{Val{0}}, differential_vars
    )
    Θ1 = 1 - Θ
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, idxs)
        if m === true
            @views @.. Θ1 * y₀[idxs] + Θ * y₁[idxs] +
                Θ * (Θ - 1) * ((1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs])
        elseif m === false
            @views @.. Θ1 * y₀[idxs] + Θ * y₁[idxs]
        else
            @views @.. Θ1 * y₀[idxs] + Θ * y₁[idxs] +
                m * Θ * (Θ - 1) * ((1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs])
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. Θ1 * y₀[idxs] + Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
    elseif cache.interp_order == 4
        @views @.. Θ1 * y₀[idxs] + Θ * (
            y₁[idxs] + Θ1 * (
                k[1][idxs] +
                    Θ * (k[2][idxs] + Θ * (k[3][idxs] + Θ * k[4][idxs]))
            )
        )
    else
        @views @.. Θ1 * y₀[idxs] +
            Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
    end
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            RosenbrockCombinedConstantCache, RosenbrockCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    Θ1 = 1 - Θ
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, nothing)
        if m === true
            @.. out = Θ1 * y₀ + Θ * y₁ +
                Θ * (Θ - 1) * ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
        elseif m === false
            @.. out = Θ1 * y₀ + Θ * y₁
        else
            @.. out = Θ1 * y₀ + Θ * y₁ +
                m * Θ * (Θ - 1) * ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. out = Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
    elseif cache.interp_order == 4
        @.. out = Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * (k[3] + Θ * k[4]))))
    else
        @.. out = Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
    end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            RosenbrockCombinedConstantCache, RosenbrockCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs, T::Type{Val{0}}, differential_vars
    )
    Θ1 = 1 - Θ
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, idxs)
        if m === true
            @views @.. out = Θ1 * y₀[idxs] + Θ * y₁[idxs] +
                Θ * (Θ - 1) * ((1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs])
        elseif m === false
            @views @.. out = Θ1 * y₀[idxs] + Θ * y₁[idxs]
        else
            @views @.. out = Θ1 * y₀[idxs] + Θ * y₁[idxs] +
                m * Θ * (Θ - 1) * ((1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs])
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. out = Θ1 * y₀[idxs] + Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
    elseif cache.interp_order == 4
        @views @.. out = Θ1 * y₀[idxs] + Θ * (
            y₁[idxs] + Θ1 * (
                k[1][idxs] +
                    Θ * (k[2][idxs] + Θ * (k[3][idxs] + Θ * k[4][idxs]))
            )
        )
    else
        @views @.. out = Θ1 * y₀[idxs] +
            Θ * (
            y₁[idxs] +
                Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs]))
        )
    end
    out
end

# First Derivative
@muladd function _ode_interpolant(
        Θ, dt,
        y₀,
        y₁,
        k,
        cache::Union{
            RosenbrockCache, RosenbrockCombinedConstantCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, nothing)
        hermite_deriv = @.. (
            k[1] + Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ * (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                    6 * y₁
            )
        ) / dt
        if m === true
            hermite_deriv
        elseif m === false
            @.. (y₁ - y₀) / dt
        else
            @.. m * hermite_deriv + (true - m) * (y₁ - y₀) / dt
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) / dt
    elseif cache.interp_order == 4
        @.. (
            k[1] + Θ * (
                -2 * k[1] + 2 * k[2] +
                    Θ * (
                    -3 * k[2] + 3 * k[3] +
                        Θ * (-4 * k[3] + 4 * k[4] - 5 * Θ * k[4])
                )
            ) -
                y₀ + y₁
        ) / dt
    else
        @.. (
            k[1] + Θ * (
                -2 * k[1] + 2 * k[2] +
                    Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])
            ) -
                y₀ + y₁
        ) / dt
    end
end
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            RosenbrockCombinedConstantCache, RosenbrockCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs, T::Type{Val{1}}, differential_vars
    )
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, idxs)
        hermite_deriv = @views @.. (
            k[1][idxs] + Θ * (
                -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] +
                    Θ * (3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] + 6 * y₀[idxs] - 6 * y₁[idxs]) +
                    6 * y₁[idxs]
            )
        ) / dt
        if m === true
            hermite_deriv
        elseif m === false
            @views @.. (y₁[idxs] - y₀[idxs]) / dt
        else
            @views @.. m * hermite_deriv + (true - m) * (y₁[idxs] - y₀[idxs]) / dt
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. (
            k[1][idxs] +
                Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] - 3 * k[2][idxs] * Θ) -
                y₀[idxs] + y₁[idxs]
        ) / dt
    elseif cache.interp_order == 4
        @views @.. (
            k[1][idxs] + Θ * (
                -2 * k[1][idxs] + 2 * k[2][idxs] +
                    Θ * (
                    -3 * k[2][idxs] + 3 * k[3][idxs] +
                        Θ * (-4 * k[3][idxs] + 4 * k[4][idxs] - 5 * Θ * k[4][idxs])
                )
            ) -
                y₀[idxs] + y₁[idxs]
        ) / dt
    else
        @views @.. (
            k[1][idxs] +
                Θ * (
                -2 * k[1][idxs] + 2 * k[2][idxs] +
                    Θ * (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])
            ) -
                y₀[idxs] + y₁[idxs]
        ) / dt
    end
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            RosenbrockCombinedConstantCache, RosenbrockCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, nothing)
        if m === true
            @.. out = (
                k[1] + Θ * (
                    -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                        Θ * (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                        6 * y₁
                )
            ) / dt
        elseif m === false
            @.. out = (y₁ - y₀) / dt
        else
            @.. out = m * (
                k[1] + Θ * (
                    -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                        Θ * (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                        6 * y₁
                )
            ) / dt + (true - m) * (y₁ - y₀) / dt
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. out = (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) / dt
    elseif cache.interp_order == 4
        @.. out = (
            k[1] + Θ * (
                -2 * k[1] + 2 * k[2] +
                    Θ * (
                    -3 * k[2] + 3 * k[3] +
                        Θ * (-4 * k[3] + 4 * k[4] - 5 * Θ * k[4])
                )
            ) -
                y₀ + y₁
        ) / dt
    else
        @.. out = (
            k[1] +
                Θ * (
                -2 * k[1] + 2 * k[2] +
                    Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])
            ) - y₀ + y₁
        ) / dt
    end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{
            RosenbrockCombinedConstantCache, RosenbrockCache,

            HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
        },
        idxs, T::Type{Val{1}}, differential_vars
    )
    if hasproperty(cache, :interp_order) && cache.interp_order == -1
        m = _dae_hermite_mask(differential_vars, idxs)
        if m === true
            @views @.. out = (
                k[1][idxs] + Θ * (
                    -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] +
                        Θ * (3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] + 6 * y₀[idxs] - 6 * y₁[idxs]) +
                        6 * y₁[idxs]
                )
            ) / dt
        elseif m === false
            @views @.. out = (y₁[idxs] - y₀[idxs]) / dt
        else
            @views @.. out = m * (
                k[1][idxs] + Θ * (
                    -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] +
                        Θ * (3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] + 6 * y₀[idxs] - 6 * y₁[idxs]) +
                        6 * y₁[idxs]
                )
            ) / dt + (true - m) * (y₁[idxs] - y₀[idxs]) / dt
        end
    elseif !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. out = (
            k[1][idxs] +
                Θ *
                (
                -2 * k[1][idxs] + 2 * k[2][idxs] -
                    3 * k[2][idxs] * Θ
            ) -
                y₀[idxs] + y₁[idxs]
        ) / dt
    elseif cache.interp_order == 4
        @views @.. out = (
            k[1][idxs] + Θ * (
                -2 * k[1][idxs] + 2 * k[2][idxs] +
                    Θ * (
                    -3 * k[2][idxs] + 3 * k[3][idxs] +
                        Θ * (-4 * k[3][idxs] + 4 * k[4][idxs] - 5 * Θ * k[4][idxs])
                )
            ) -
                y₀[idxs] + y₁[idxs]
        ) / dt
    else
        @views @.. broadcast = false out = (
            k[1][idxs] +
                Θ * (
                -2 * k[1][idxs] + 2 * k[2][idxs] +
                    Θ *
                    (
                    -3 * k[2][idxs] + 3 * k[3][idxs] -
                        4 * Θ * k[3][idxs]
                )
            ) -
                y₀[idxs] + y₁[idxs]
        ) / dt
    end
    out
end

# Second Derivative
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::RosenbrockCombinedConstantCache,
        idxs::Nothing, T::Type{Val{2}}, differential_vars
    )
    if cache.interp_order == 4
        @inbounds (
            -2 * k[1] + 2 * k[2] + Θ * (
                -6 * k[2] + 6 * k[3] +
                    Θ * (-12 * k[3] + 12 * k[4] - 20 * Θ * k[4])
            )
        ) / dt^2
    else
        @inbounds (-2 * k[1] + 2 * k[2] + Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])) / dt^2
    end
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::RosenbrockCache, idxs::Nothing,
        T::Type{Val{2}}, differential_vars
    )
    if cache.interp_order == 4
        @inbounds @.. broadcast = false (
            -2 * k[1] + 2 * k[2] + Θ * (
                -6 * k[2] + 6 * k[3] +
                    Θ * (-12 * k[3] + 12 * k[4] - 20 * Θ * k[4])
            )
        ) / dt^2
    else
        @inbounds @.. broadcast = false (
            -2 * k[1] + 2 * k[2] +
                Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])
        ) / dt^2
    end
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache, HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache},
        idxs, T::Type{Val{2}}, differential_vars
    )
    if cache.interp_order == 4
        @.. broadcast = false (
            -2 * k[1][idxs] + 2 * k[2][idxs] + Θ * (
                -6 * k[2][idxs] + 6 * k[3][idxs] +
                    Θ * (-12 * k[3][idxs] + 12 * k[4][idxs] - 20 * Θ * k[4][idxs])
            )
        ) / dt^2
    else
        @.. broadcast = false (
            -2 * k[1][idxs] + 2 * k[2][idxs] +
                Θ * (-6 * k[2][idxs] + 6 * k[3][idxs] - 12 * Θ * k[3][idxs])
        ) / dt^2
    end
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache, HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars
    )
    if cache.interp_order == 4
        @.. broadcast = false out = (
            -2 * k[1] + 2 * k[2] + Θ * (
                -6 * k[2] + 6 * k[3] +
                    Θ * (-12 * k[3] + 12 * k[4] - 20 * Θ * k[4])
            )
        ) / dt^2
    else
        @.. broadcast = false out = (
            -2 * k[1] + 2 * k[2] +
                Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])
        ) / dt^2
    end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache, HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache},
        idxs, T::Type{Val{2}}, differential_vars
    )
    if cache.interp_order == 4
        @views @.. broadcast = false out = (
            -2 * k[1][idxs] + 2 * k[2][idxs] + Θ * (
                -6 * k[2][idxs] + 6 * k[3][idxs] +
                    Θ * (-12 * k[3][idxs] + 12 * k[4][idxs] - 20 * Θ * k[4][idxs])
            )
        ) / dt^2
    else
        @views @.. broadcast = false out = (
            -2 * k[1][idxs] + 2 * k[2][idxs] +
                Θ *
                (-6 * k[2][idxs] + 6 * k[3][idxs] - 12 * Θ * k[3][idxs])
        ) /
            dt^2
    end
    out
end

# Third Derivative
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::RosenbrockCombinedConstantCache,
        idxs::Nothing, T::Type{Val{3}}, differential_vars
    )
    if cache.interp_order == 4
        @inbounds (-6 * k[2] + 6 * k[3] + Θ * (-24 * k[3] + 24 * k[4] - 60 * Θ * k[4])) / dt^3
    else
        @inbounds (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
    end
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::RosenbrockCache, idxs::Nothing,
        T::Type{Val{3}}, differential_vars
    )
    if cache.interp_order == 4
        @inbounds @.. broadcast = false (-6 * k[2] + 6 * k[3] + Θ * (-24 * k[3] + 24 * k[4] - 60 * Θ * k[4])) / dt^3
    else
        @inbounds @.. broadcast = false (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
    end
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache, HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache},
        idxs, T::Type{Val{3}}, differential_vars
    )
    if cache.interp_order == 4
        @.. broadcast = false (
            -6 * k[2][idxs] + 6 * k[3][idxs] +
                Θ * (-24 * k[3][idxs] + 24 * k[4][idxs] - 60 * Θ * k[4][idxs])
        ) / dt^3
    else
        @.. broadcast = false (-6 * k[2][idxs] + 6 * k[3][idxs] - 24 * Θ * k[3][idxs]) / dt^3
    end
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache, HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars
    )
    if cache.interp_order == 4
        @.. broadcast = false out = (-6 * k[2] + 6 * k[3] + Θ * (-24 * k[3] + 24 * k[4] - 60 * Θ * k[4])) / dt^3
    else
        @.. broadcast = false out = (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
    end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache, HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache},
        idxs, T::Type{Val{3}}, differential_vars
    )
    if cache.interp_order == 4
        @views @.. broadcast = false out = (
            -6 * k[2][idxs] + 6 * k[3][idxs] +
                Θ * (-24 * k[3][idxs] + 24 * k[4][idxs] - 60 * Θ * k[4][idxs])
        ) / dt^3
    else
        @views @.. broadcast = false out = (
            -6 * k[2][idxs] + 6 * k[3][idxs] -
                24 * Θ * k[3][idxs]
        ) / dt^3
    end
    out
end
