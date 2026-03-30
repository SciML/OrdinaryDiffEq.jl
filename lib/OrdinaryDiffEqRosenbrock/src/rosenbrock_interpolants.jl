### Fallbacks to capture
ROSENBROCKS_WITH_INTERPOLATIONS = Union{
    RosenbrockCombinedConstantCache,
    RosenbrockCache,
    HybridExplicitImplicitConstantCache, HybridExplicitImplicitCache,
}

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
        # Generic Hermite: k[1]=f₀, k[2]=f₁ (methods with no H matrix)
        @.. Θ1 * y₀ + Θ * y₁ +
            Θ * (Θ - 1) * ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
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
        @views @.. Θ1 * y₀[idxs] + Θ * y₁[idxs] +
            Θ * (Θ - 1) * ((1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs])
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
        @.. out = Θ1 * y₀ + Θ * y₁ +
            Θ * (Θ - 1) * ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
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
        @views @.. out = Θ1 * y₀[idxs] + Θ * y₁[idxs] +
            Θ * (Θ - 1) * ((1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs])
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
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
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
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
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
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
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
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
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
