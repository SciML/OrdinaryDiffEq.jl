### Fallbacks to capture

ROSENBROCKS_WITH_INTERPOLATIONS = Union{Rosenbrock23ConstantCache, Rosenbrock23Cache,
    Rosenbrock32ConstantCache, Rosenbrock32Cache,
    Rodas23WConstantCache, Rodas3PConstantCache,
    Rodas23WCache, Rodas3PCache,
    Rodas4ConstantCache, Rosenbrock5ConstantCache,
    Rodas4Cache, Rosenbrock5Cache}

function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::ROSENBROCKS_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::ROSENBROCKS_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
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

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars)
    @rosenbrock2332pre0
    @inbounds y₀ + dt * (c1 * k[1] + c2 * k[2])
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23Cache, Rosenbrock32Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @rosenbrock2332pre0
    @inbounds @.. broadcast=false y₀+dt * (c1 * k[1] + c2 * k[2])
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache, Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs, T::Type{Val{0}}, differential_vars)
    @rosenbrock2332pre0
    @.. broadcast=false y₀[idxs]+dt * (c1 * k[1][idxs] + c2 * k[2][idxs])
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs::Nothing, T::Type{Val{0}}, differential_vars)
    @rosenbrock2332pre0
    @inbounds @.. broadcast=false out=y₀ + dt * (c1 * k[1] + c2 * k[2])
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs, T::Type{Val{0}}, differential_vars)
    @rosenbrock2332pre0
    @views @.. broadcast=false out=y₀[idxs] + dt * (c1 * k[1][idxs] + c2 * k[2][idxs])
    out
end

# First Derivative of the dense output
@def rosenbrock2332pre1 begin
    @rosenbrock2332unpack
    c1diff = (1 - 2 * Θ) / (1 - 2 * d)
    c2diff = (2 * Θ - 2 * d) / (1 - 2 * d)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache, Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @rosenbrock2332pre1
    @.. broadcast=false c1diff * k[1]+c2diff * k[2]
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache, Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs, T::Type{Val{1}}, differential_vars)
    @rosenbrock2332pre1
    @.. broadcast=false c1diff * k[1][idxs]+c2diff * k[2][idxs]
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @rosenbrock2332pre1
    @.. broadcast=false out=c1diff * k[1] + c2diff * k[2]
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock23Cache,
            Rosenbrock32ConstantCache, Rosenbrock32Cache
        }, idxs, T::Type{Val{1}}, differential_vars)
    @rosenbrock2332pre1
    @views @.. broadcast=false out=c1diff * k[1][idxs] + c2diff * k[2][idxs]
    out
end

"""
From MATLAB ODE Suite by Shampine
"""
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas23WConstantCache, Rodas3PConstantCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @inbounds Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{Rodas4Cache, Rodas23WCache, Rodas3PCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @inbounds @.. broadcast=false Θ1 * y₀+Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas4Cache, Rodas23WConstantCache,
            Rodas23WCache, Rodas3PConstantCache, Rodas3PCache},
        idxs, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @.. broadcast=false Θ1 * y₀[idxs]+Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas4Cache, Rodas23WConstantCache,
            Rodas23WCache, Rodas3PConstantCache, Rodas3PCache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @.. broadcast=false out=Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas4Cache, Rodas23WConstantCache,
            Rodas23WCache, Rodas3PConstantCache, Rodas3PCache},
        idxs, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @views @.. broadcast=false out=Θ1 * y₀[idxs] +
                                   Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
    out
end

# First Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas23WConstantCache, Rodas3PConstantCache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @inbounds (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) / dt
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Union{Rodas4Cache, Rodas23WCache, Rodas3PCache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @inbounds @.. broadcast=false (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ +
                                   y₁)/dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas4Cache, Rodas23WConstantCache,
            Rodas23WCache, Rodas3PConstantCache, Rodas3PCache},
        idxs, T::Type{Val{1}}, differential_vars)
    @.. broadcast=false (k[1][idxs] +
                         Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] - 3 * k[2][idxs] * Θ) -
                         y₀[idxs] + y₁[idxs])/dt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas4Cache, Rodas23WConstantCache,
            Rodas23WCache, Rodas3PConstantCache, Rodas3PCache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @.. broadcast=false out=(k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) /
                            dt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rodas4ConstantCache, Rodas4Cache, Rodas23WConstantCache,
            Rodas23WCache, Rodas3PConstantCache, Rodas3PCache},
        idxs, T::Type{Val{1}}, differential_vars)
    @views @.. broadcast=false out=(k[1][idxs] +
                                    Θ *
                                    (-2 * k[1][idxs] + 2 * k[2][idxs] -
                                     3 * k[2][idxs] * Θ) -
                                    y₀[idxs] + y₁[idxs]) / dt
    out
end

#-

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5ConstantCache,
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @inbounds Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @inbounds @.. broadcast=false Θ1 * y₀+Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @.. broadcast=false Θ1 *
                        y₀[idxs]+Θ * (y₁[idxs] +
                                  Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @.. broadcast=false out=Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{0}}, differential_vars)
    Θ1 = 1 - Θ
    @views @.. broadcast=false out=Θ1 * y₀[idxs] +
                                   Θ * (y₁[idxs] +
                                    Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
    out
end

# First Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5ConstantCache,
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @inbounds (k[1] +
               Θ * (-2 * k[1] + 2 * k[2] + Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ +
               y₁) / dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache, idxs::Nothing,
        T::Type{Val{1}}, differential_vars)
    @inbounds @.. broadcast=false (k[1] +
                                   Θ * (-2 * k[1] + 2 * k[2] +
                                    Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ + y₁)/dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{1}}, differential_vars)
    @.. broadcast=false (k[1][idxs] +
                         Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] +
                          Θ * (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])) -
                         y₀[idxs] + y₁[idxs])/dt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars)
    @.. broadcast=false out=(k[1] +
                             Θ * (-2 * k[1] + 2 * k[2] +
                              Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ + y₁) / dt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{1}}, differential_vars)
    @views @.. broadcast=false out=(k[1][idxs] +
                                    Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] +
                                     Θ *
                                     (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])) -
                                    y₀[idxs] + y₁[idxs]) / dt
    out
end

# Second Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5ConstantCache,
        idxs::Nothing, T::Type{Val{2}}, differential_vars)
    @inbounds (-2 * k[1] + 2 * k[2] + Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])) / dt^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache, idxs::Nothing,
        T::Type{Val{2}}, differential_vars)
    @inbounds @.. broadcast=false (-2 * k[1] + 2 * k[2] +
                                   Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3]))/dt^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{2}}, differential_vars)
    @.. broadcast=false (-2 * k[1][idxs] + 2 * k[2][idxs] +
                         Θ * (-6 * k[2][idxs] + 6 * k[3][idxs] - 12 * Θ * k[3][idxs]))/dt^2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars)
    @.. broadcast=false out=(-2 * k[1] + 2 * k[2] +
                             Θ * (-6 * k[2] + 6 * k[3] - 12 * Θ * k[3])) / dt^2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{2}}, differential_vars)
    @views @.. broadcast=false out=(-2 * k[1][idxs] + 2 * k[2][idxs] +
                                    Θ *
                                    (-6 * k[2][idxs] + 6 * k[3][idxs] - 12 * Θ * k[3][idxs])) /
                                   dt^2
    out
end

# Third Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5ConstantCache,
        idxs::Nothing, T::Type{Val{3}}, differential_vars)
    @inbounds (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache, idxs::Nothing,
        T::Type{Val{3}}, differential_vars)
    @inbounds @.. broadcast=false (-6 * k[2] + 6 * k[3] - 24 * Θ * k[3])/dt^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{3}}, differential_vars)
    @.. broadcast=false (-6 * k[2][idxs] + 6 * k[3][idxs] - 24 * Θ * k[3][idxs])/dt^3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars)
    @.. broadcast=false out=(-6 * k[2] + 6 * k[3] - 24 * Θ * k[3]) / dt^3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
        idxs, T::Type{Val{3}}, differential_vars)
    @views @.. broadcast=false out=(-6 * k[2][idxs] + 6 * k[3][idxs] -
                                    24 * Θ * k[3][idxs]) / dt^3
    out
end
