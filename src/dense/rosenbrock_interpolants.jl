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
    T::Type{Val{0}}, dv=nothing)
    @rosenbrock2332pre0
    @inbounds y₀ + dt * (c1 * k[1] + c2 * k[2])
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23Cache, Rosenbrock32Cache},
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    @rosenbrock2332pre0
    @inbounds @.. broadcast=false y₀+dt * (c1 * k[1] + c2 * k[2])
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23ConstantCache, Rosenbrock23Cache,
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs, T::Type{Val{0}}, dv=nothing)
    @rosenbrock2332pre0
    @.. broadcast=false y₀[idxs]+dt * (c1 * k[1][idxs] + c2 * k[2][idxs])
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23ConstantCache,
        Rosenbrock23Cache,
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    @rosenbrock2332pre0
    @inbounds @.. broadcast=false out=y₀ + dt * (c1 * k[1] + c2 * k[2])
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k, cache::Rosenbrock23Cache{<:Array},
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    @rosenbrock2332pre0
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = y₀[i] + dt * (c1 * k[1][i] + c2 * k[2][i])
    end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23ConstantCache,
        Rosenbrock23Cache,
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs, T::Type{Val{0}}, dv=nothing)
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
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs::Nothing, T::Type{Val{1}}, dv=nothing)
    @rosenbrock2332pre1
    @.. broadcast=false c1diff * k[1]+c2diff * k[2]
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23ConstantCache, Rosenbrock23Cache,
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs, T::Type{Val{1}}, dv=nothing)
    @rosenbrock2332pre1
    @.. broadcast=false c1diff * k[1][idxs]+c2diff * k[2][idxs]
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23ConstantCache,
        Rosenbrock23Cache,
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs::Nothing, T::Type{Val{1}}, dv=nothing)
    @rosenbrock2332pre1
    @.. broadcast=false out=c1diff * k[1] + c2diff * k[2]
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock23ConstantCache,
        Rosenbrock23Cache,
        Rosenbrock32ConstantCache, Rosenbrock32Cache,
    }, idxs, T::Type{Val{1}}, dv=nothing)
    @rosenbrock2332pre1
    @views @.. broadcast=false out=c1diff * k[1][idxs] + c2diff * k[2][idxs]
    out
end

"""
From MATLAB ODE Suite by Shampine
"""
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rodas4ConstantCache,
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @inbounds Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rodas4Cache, idxs::Nothing,
    T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @inbounds @.. broadcast=false Θ1 * y₀+Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rodas4ConstantCache, Rodas4Cache}, idxs,
    T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @.. broadcast=false Θ1 * y₀[idxs]+Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rodas4ConstantCache, Rodas4Cache},
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @.. broadcast=false out=Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * k[2]))
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k, cache::Rodas4Cache{<:Array},
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = Θ1 * y₀[i] + Θ * (y₁[i] + Θ1 * (k[1][i] + Θ * k[2][i]))
    end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rodas4ConstantCache, Rodas4Cache}, idxs,
    T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @views @.. broadcast=false out=Θ1 * y₀[idxs] +
                                   Θ * (y₁[idxs] + Θ1 * (k[1][idxs] + Θ * k[2][idxs]))
    out
end

# First Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rodas4ConstantCache,
    idxs::Nothing, T::Type{Val{1}}, dv=nothing)
    @inbounds (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) / dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rodas4Cache, idxs::Nothing,
    T::Type{Val{1}}, dv=nothing)
    @inbounds @.. broadcast=false (k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ +
                                   y₁)/dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rodas4ConstantCache, Rodas4Cache}, idxs,
    T::Type{Val{1}}, dv=nothing)
    @.. broadcast=false (k[1][idxs] +
                         Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] - 3 * k[2][idxs] * Θ) -
                         y₀[idxs] + y₁[idxs])/dt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rodas4ConstantCache, Rodas4Cache},
    idxs::Nothing, T::Type{Val{1}}, dv=nothing)
    @.. broadcast=false out=(k[1] + Θ * (-2 * k[1] + 2 * k[2] - 3 * k[2] * Θ) - y₀ + y₁) /
                            dt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rodas4ConstantCache, Rodas4Cache}, idxs,
    T::Type{Val{1}}, dv=nothing)
    @views @.. broadcast=false out=(k[1][idxs] +
                                    Θ *
                                    (-2 * k[1][idxs] + 2 * k[2][idxs] -
                                     3 * k[2][idxs] * Θ) -
                                    y₀[idxs] + y₁[idxs]) / dt
    out
end

#-

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5ConstantCache,
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @inbounds Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache, idxs::Nothing,
    T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @inbounds @.. broadcast=false Θ1 * y₀+Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
    idxs, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @.. broadcast=false Θ1 *
                        y₀[idxs]+Θ * (y₁[idxs] +
                                  Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @.. broadcast=false out=Θ1 * y₀ + Θ * (y₁ + Θ1 * (k[1] + Θ * (k[2] + Θ * k[3])))
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache{<:Array},
    idxs::Nothing, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = Θ1 * y₀[i] + Θ * (y₁[i] + Θ1 * (k[1][i] + Θ * (k[2][i] + Θ * k[3][i])))
    end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
    idxs, T::Type{Val{0}}, dv=nothing)
    Θ1 = 1 - Θ
    @views @.. broadcast=false out=Θ1 * y₀[idxs] +
                                   Θ * (y₁[idxs] +
                                    Θ1 * (k[1][idxs] + Θ * (k[2][idxs] + Θ * k[3][idxs])))
    out
end

# First Derivative
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5ConstantCache,
    idxs::Nothing, T::Type{Val{1}}, dv=nothing)
    @inbounds (k[1] +
               Θ * (-2 * k[1] + 2 * k[2] + Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ +
               y₁) / dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Rosenbrock5Cache, idxs::Nothing,
    T::Type{Val{1}}, dv=nothing)
    @inbounds @.. broadcast=false (k[1] +
                                   Θ * (-2 * k[1] + 2 * k[2] +
                                    Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ + y₁)/dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
    idxs, T::Type{Val{1}}, dv=nothing)
    @.. broadcast=false (k[1][idxs] +
                         Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] +
                          Θ * (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])) -
                         y₀[idxs] + y₁[idxs])/dt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
    idxs::Nothing, T::Type{Val{1}}, dv=nothing)
    @.. broadcast=false out=(k[1] +
                             Θ * (-2 * k[1] + 2 * k[2] +
                              Θ * (-3 * k[2] + 3 * k[3] - 4 * Θ * k[3])) - y₀ + y₁) / dt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{Rosenbrock5ConstantCache, Rosenbrock5Cache},
    idxs, T::Type{Val{1}}, dv=nothing)
    @views @.. broadcast=false out=(k[1][idxs] +
                                    Θ * (-2 * k[1][idxs] + 2 * k[2][idxs] +
                                     Θ *
                                     (-3 * k[2][idxs] + 3 * k[3][idxs] - 4 * Θ * k[3][idxs])) -
                                    y₀[idxs] + y₁[idxs]) / dt
    out
end
#-
