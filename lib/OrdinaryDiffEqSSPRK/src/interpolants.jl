const NEGZERO = Float16(-0.0f0)

@def ssprkpre0 begin
    c00 = @evalpoly(Θ, 1, NEGZERO, -1)
    c10 = Θ^2
    b10dt = Θ * @evalpoly(Θ, 1, -1) * dt
end

@def ssprkpre1 begin
    b10diff = @evalpoly(Θ, 1, -2)
    c10diffinvdt = 2Θ * inv(dt) # = -c00diff * inv(dt)
end

@def ssprkpre2 begin
    invdt = inv(dt)
    b10diff2invdt = -2 * invdt
    c10diff2invdt2 = 2 * invdt^2 # = -c00diff2 * inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
@ssprkpre0
@inbounds @.. broadcast=false y₀*c00+y₁*c10+k[1]*b10dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache}, idxs,
    T::Type{Val{0}}, differential_vars::Nothing)
@ssprkpre0
@views @.. broadcast=false y₀[idxs]*c00+y₁[idxs]*c10+k[1][idxs]*b10dt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache},
    idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
@ssprkpre0
@inbounds @.. broadcast=false out=y₀ * c00 + y₁ * c10 + k[1] * b10dt
out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache}, idxs,
    T::Type{Val{0}}, differential_vars::Nothing)
@ssprkpre0
@views @.. broadcast=false out=y₀[idxs] * c00 + y₁[idxs] * c10 + k[1][idxs] * b10dt
out
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache},
    idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
@ssprkpre1
@inbounds @.. broadcast=false (y₁ - y₀) * c10diffinvdt+k[1] * b10diff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache}, idxs,
    T::Type{Val{1}}, differential_vars::Nothing)
@ssprkpre1
@views @.. broadcast=false (y₁[idxs] - y₀[idxs]) * c10diffinvdt+k[1][idxs] * b10diff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache},
    idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
@ssprkpre1
@inbounds @.. broadcast=false out=(y₁ - y₀) * c10diffinvdt + k[1] * b10diff
out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
            SSPRK33ConstantCache, SSPRK33Cache,
            SSPRK43ConstantCache, SSPRK43Cache,
            SSPRK432ConstantCache, SSPRK432Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    @ssprkpre1
    @views @.. broadcast=false out=(y₁[idxs] - y₀[idxs]) * c10diffinvdt +
                                   k[1][idxs] * b10diff
    out
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
            SSPRK33ConstantCache, SSPRK33Cache,
            SSPRK43ConstantCache, SSPRK43Cache,
            SSPRK432ConstantCache, SSPRK432Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @ssprkpre2
    @inbounds @.. broadcast=false (y₁ - y₀) * c10diff2invdt2+k[1] * b10diff2invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache}, idxs,
    T::Type{Val{2}}, differential_vars::Nothing)
@ssprkpre2
@views @.. broadcast=false (y₁[idxs] - y₀[idxs]) *
                           c10diff2invdt2+k[1][idxs] * b10diff2invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
    cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
        SSPRK33ConstantCache, SSPRK33Cache,
        SSPRK43ConstantCache, SSPRK43Cache,
        SSPRK432ConstantCache, SSPRK432Cache},
    idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
@ssprkpre2
@inbounds @.. broadcast=false out=(y₁ - y₀) * c10diff2invdt2 + k[1] * b10diff2invdt
out
end
