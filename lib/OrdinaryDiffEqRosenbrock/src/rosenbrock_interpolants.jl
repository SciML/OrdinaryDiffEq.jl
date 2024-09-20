"""
From MATLAB ODE Suite by Shampine
"""

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. Θ1 * y₀+Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
    else
        @.. Θ1 * y₀+Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
    end
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. Θ1 * y₀[idxs]+Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
    else
        @views @.. Θ1 * y₀[idxs]+Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
    end
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. out=Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
    else
        @.. out=Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
    end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. out=Θ1 * y₀[idxs] + Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
    else
        @views @.. out=Θ1 * y₀[idxs]+Θ * (y₁[idxs] +
            Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
    end
    out
end

# First Derivative
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{RosenbrockCache, RosenbrockCombinedConstantCache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁)/dt
    else
        @.. (k[1] +  Θ * (-2 * k[1] + 2 * k[2] +
             Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) -
             y₀ + y₁)/dt
    end
end
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{1}}, differential_vars)
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. (k[1][idxs] + Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] - 3 * k[2][idxs] * Θ) -
                    y₀[idxs] + y₁[idxs])/dt
    else
       @views  @.. (k[1][idxs] + Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] +
             Θ * (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])) - y₀[idxs] + y₁[idxs])/dt
    end
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @.. out=(k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) / dt
    else
        @.. out=(k[1] + Θ * (-2 * k[1] + 2 * k[2] +
                Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ + y₁) / dt
    end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{1}}, differential_vars)
    if !hasproperty(cache, :interp_order) || cache.interp_order == 2
        @views @.. out=(k[1][idxs] +
                        Θ *
                        (-2 * k[1][idxs] + 2 * k[2][idxs] -
                         3 * k[2][idxs] * Θ) -
                        y₀[idxs] + y₁[idxs]) / dt
    else
        @views @.. broadcast=false out=(k[1][idxs] +
                                        Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] +
                                         Θ *
                                         (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])) -
                                        y₀[idxs] + y₁[idxs]) / dt
    end
    out
end

# Second Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::RosenbrockCombinedConstantCache,
        idxs::Nothing, T::Type{Val{2}}, differential_vars)
    @inbounds (-2 * k[1] + 2 * k[2] + Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])) / dt^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::RosenbrockCache, idxs::Nothing,
        T::Type{Val{2}}, differential_vars)
    @inbounds @.. broadcast=false (-2 * k[1] + 2 * k[2] +
                                   Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3]))/dt^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{2}}, differential_vars)
    @.. broadcast=false (-2 * k[1][idxs] + 2 * k[2][idxs] +
                         Θ * (-6 * k[2][idxs] + 6 * k[3][idxs] - 12 * Θ * k[3][idxs]))/dt^2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars)
    @.. broadcast=false out=(-2 * k[1] + 2 * k[2] +
                             Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])) / dt^2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{2}}, differential_vars)
    @views @.. broadcast=false out=(-2 * k[1][idxs] + 2 * k[2][idxs] +
                                    Θ *
                                    (-6 * k[2][idxs] + 6 * k[3][idxs] - 12 * Θ * k[3][idxs])) /
                                   dt^2
    out
end

# Third Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::RosenbrockCombinedConstantCache,
        idxs::Nothing, T::Type{Val{3}}, differential_vars)
    @inbounds (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::RosenbrockCache, idxs::Nothing,
        T::Type{Val{3}}, differential_vars)
    @inbounds @.. broadcast=false (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3])/dt^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{3}}, differential_vars)
    @.. broadcast=false (-6 * k[2][idxs] + 6 * k[3][idxs] - 24 * Θ * k[3][idxs])/dt^3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars)
    @.. broadcast=false out=(-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{RosenbrockCombinedConstantCache, RosenbrockCache},
        idxs, T::Type{Val{3}}, differential_vars)
    @views @.. broadcast=false out=(-6 * k[2][idxs] + 6 * k[3][idxs] -
                                    24 * Θ * k[3][idxs]) / dt^3
    out
end
