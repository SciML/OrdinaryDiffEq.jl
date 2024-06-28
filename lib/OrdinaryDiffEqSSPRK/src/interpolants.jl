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

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
            SSPRK33ConstantCache, SSPRK33Cache,
            SSPRK43ConstantCache, SSPRK43Cache,
            SSPRK432ConstantCache, SSPRK432Cache}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing)
    @ssprkpre2
    @views @.. broadcast=false out=(y₁[idxs] - y₀[idxs]) * c10diff2invdt2 +
                                   k[1][idxs] * b10diff2invdt
    out
end