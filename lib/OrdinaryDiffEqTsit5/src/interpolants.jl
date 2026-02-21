RK_WITH_SPECIAL_INTERPOLATIONS = Union{Tsit5ConstantCache, Tsit5CacheType}

function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::RK_WITH_SPECIAL_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::RK_WITH_SPECIAL_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

"""
Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
@def tsit5unpack begin
    var"#T#" = constvalue(recursive_unitless_bottom_eltype(y₁))
    r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62, r63, r64, r72, r73, r74 = Tsit5Interp(var"#T#")
end

@def tsit5pre0 begin
    @tsit5unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, r11, r12, r13, r14)
    b2Θ = Θ² * @evalpoly(Θ, r22, r23, r24)
    b3Θ = Θ² * @evalpoly(Θ, r32, r33, r34)
    b4Θ = Θ² * @evalpoly(Θ, r42, r43, r44)
    b5Θ = Θ² * @evalpoly(Θ, r52, r53, r54)
    b6Θ = Θ² * @evalpoly(Θ, r62, r63, r64)
    b7Θ = Θ² * @evalpoly(Θ, r72, r73, r74)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Tsit5ConstantCache,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[2]*b2Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ)
    return @inbounds y₀ +
        dt * (
        k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ +
            k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Tsit5CacheType, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    return @inbounds @.. broadcast = false y₀ + dt * (
        k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ +
            k[4] * b4Θ +
            k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    return y₀[idxs] +
        dt * (
        k[1][idxs] * b1Θ + k[2][idxs] * b2Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ
    )
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    @inbounds @.. broadcast = false out = y₀ +
        dt *
        (
        k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ +
            k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ
    )
    out
end

@muladd function _ode_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = y₀[i] +
            dt * (
            k[1][i] * b1Θ + k[2][i] * b2Θ + k[3][i] * b3Θ + k[4][i] * b4Θ +
                k[5][i] * b5Θ + k[6][i] * b6Θ + k[7][i] * b7Θ
        )
    end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    @views @.. broadcast = false out = y₀[idxs] +
        dt *
        (
        k[1][idxs] * b1Θ + k[2][idxs] * b2Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ +
            k[7][idxs] * b7Θ
    )
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[2][i]*b2Θ + k[3][i]*b3Θ + k[4][i]*b4Θ + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ)
    #end
    out
end

@muladd function _ode_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @tsit5pre0
    @inbounds for (j, i) in enumerate(idxs)
        out[j] = y₀[i] +
            dt * (
            k[1][i] * b1Θ + k[2][i] * b2Θ + k[3][i] * b3Θ + k[4][i] * b4Θ +
                k[5][i] * b5Θ + k[6][i] * b6Θ + k[7][i] * b7Θ
        )
    end
    out
end

@def tsit5pre1 begin
    @tsit5unpack
    b1Θdiff = @evalpoly(Θ, r11, 2 * r12, 3 * r13, 4 * r14)
    b2Θdiff = Θ * @evalpoly(Θ, 2 * r22, 3 * r23, 4 * r24)
    b3Θdiff = Θ * @evalpoly(Θ, 2 * r32, 3 * r33, 4 * r34)
    b4Θdiff = Θ * @evalpoly(Θ, 2 * r42, 3 * r43, 4 * r44)
    b5Θdiff = Θ * @evalpoly(Θ, 2 * r52, 3 * r53, 4 * r54)
    b6Θdiff = Θ * @evalpoly(Θ, 2 * r62, 3 * r63, 4 * r64)
    b7Θdiff = Θ * @evalpoly(Θ, 2 * r72, 3 * r73, 4 * r74)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Tsit5ConstantCache,
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @tsit5pre1
    # return @.. broadcast=false k[1]*b1Θdiff + k[2]*b2Θdiff + k[3]*b3Θdiff + k[4]*b4Θdiff + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff
    return @inbounds k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
        k[5] * b5Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Tsit5CacheType, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @tsit5pre1
    return @inbounds @.. broadcast = false k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff +
        k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @tsit5pre1
    # return @.. broadcast=false k[1][idxs]*b1Θdiff + k[2][idxs]*b2Θdiff + k[3][idxs]*b3Θdiff + k[4][idxs]*b4Θdiff + k[5][idxs]*b5Θdiff + k[6][idxs]*b6Θdiff + k[7][idxs]*b7Θdiff
    return k[1][idxs] * b1Θdiff + k[2][idxs] * b2Θdiff + k[3][idxs] * b3Θdiff +
        k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
        k[7][idxs] * b7Θdiff
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @tsit5pre1
    @inbounds @.. broadcast = false out = k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff +
        k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
        k[7] * b7Θdiff
    #@inbounds for i in eachindex(out)
    #  out[i] = k[1][i]*b1Θdiff + k[2][i]*b2Θdiff + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff + k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff
    #end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @tsit5pre1
    @views @.. broadcast = false out = k[1][idxs] * b1Θdiff + k[2][idxs] * b2Θdiff +
        k[3][idxs] * b3Θdiff + k[4][idxs] * b4Θdiff +
        k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
        k[7][idxs] * b7Θdiff
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = k[1][i]*b1Θdiff + k[2][i]*b2Θdiff + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff + k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff
    #end
    out
end

@def tsit5pre2 begin
    @tsit5unpack
    b1Θdiff2 = @evalpoly(Θ, 2 * r12, 6 * r13, 12 * r14)
    b2Θdiff2 = @evalpoly(Θ, 2 * r22, 6 * r23, 12 * r24)
    b3Θdiff2 = @evalpoly(Θ, 2 * r32, 6 * r33, 12 * r34)
    b4Θdiff2 = @evalpoly(Θ, 2 * r42, 6 * r43, 12 * r44)
    b5Θdiff2 = @evalpoly(Θ, 2 * r52, 6 * r53, 12 * r54)
    b6Θdiff2 = @evalpoly(Θ, 2 * r62, 6 * r63, 12 * r64)
    b7Θdiff2 = @evalpoly(Θ, 2 * r72, 6 * r73, 12 * r74)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing
    )
    @tsit5pre2
    # return @.. broadcast=false k[1]*b1Θdiff2 + k[2]*b2Θdiff2 + k[3]*b3Θdiff2 + k[4]*b4Θdiff2 + k[5]*b5Θdiff2 + k[6]*b6Θdiff2 + k[7]*b7Θdiff2
    return @inbounds (
        k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
            k[4] * b4Θdiff2 +
            k[5] * b5Θdiff2 + k[6] * b6Θdiff2 + k[7] * b7Θdiff2
    ) * invdt
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing
    )
    @tsit5pre2
    # return @.. broadcast=false k[1][idxs]*b1Θdiff2 + k[2][idxs]*b2Θdiff2 + k[3][idxs]*b3Θdiff2 + k[4][idxs]*b4Θdiff2 + k[5][idxs]*b5Θdiff2 + k[6][idxs]*b6Θdiff2 + k[7][idxs]*b7Θdiff2
    return (
        k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 + k[3][idxs] * b3Θdiff2 +
            k[4][idxs] * b4Θdiff2 + k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2 +
            k[7][idxs] * b7Θdiff2
    ) * invdt
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing
    )
    @tsit5pre2
    @inbounds @.. broadcast = false out = (
        k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
            k[4] * b4Θdiff2 + k[5] * b5Θdiff2 + k[6] * b6Θdiff2 +
            k[7] * b7Θdiff2
    ) * invdt
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff2 + k[2][i]*b2Θdiff2 + k[3][i]*b3Θdiff2 + k[4][i]*b4Θdiff2 + k[5][i]*b5Θdiff2 + k[6][i]*b6Θdiff2 + k[7][i]*b7Θdiff2)*invdt
    #end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing
    )
    @tsit5pre2
    @views @.. broadcast = false out = (
        k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 +
            k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2 +
            k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2 +
            k[7][idxs] * b7Θdiff2
    ) * invdt
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff2 + k[2][i]*b2Θdiff2 + k[3][i]*b3Θdiff2 + k[4][i]*b4Θdiff2 + k[5][i]*b5Θdiff2 + k[6][i]*b6Θdiff2 + k[7][i]*b7Θdiff2)*invdt
    #end
    out
end

@def tsit5pre3 begin
    @tsit5unpack
    b1Θdiff3 = @evalpoly(Θ, 6 * r13, 24 * r14)
    b2Θdiff3 = @evalpoly(Θ, 6 * r23, 24 * r24)
    b3Θdiff3 = @evalpoly(Θ, 6 * r33, 24 * r34)
    b4Θdiff3 = @evalpoly(Θ, 6 * r43, 24 * r44)
    b5Θdiff3 = @evalpoly(Θ, 6 * r53, 24 * r54)
    b6Θdiff3 = @evalpoly(Θ, 6 * r63, 24 * r64)
    b7Θdiff3 = @evalpoly(Θ, 6 * r73, 24 * r74)
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing
    )
    @tsit5pre3
    # return @.. broadcast=false k[1]*b1Θdiff3 + k[2]*b2Θdiff3 + k[3]*b3Θdiff3 + k[4]*b4Θdiff3 + k[5]*b5Θdiff3 + k[6]*b6Θdiff3 + k[7]*b7Θdiff3
    return @inbounds (
        k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
            k[4] * b4Θdiff3 +
            k[5] * b5Θdiff3 + k[6] * b6Θdiff3 + k[7] * b7Θdiff3
    ) * invdt2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{3}}, differential_vars::Nothing
    )
    @tsit5pre3
    # return @.. broadcast=false k[1][idxs]*b1Θdiff3 + k[2][idxs]*b2Θdiff3 + k[3][idxs]*b3Θdiff3 + k[4][idxs]*b4Θdiff3 + k[5][idxs]*b5Θdiff3 + k[6][idxs]*b6Θdiff3 + k[7][idxs]*b7Θdiff3
    return (
        k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 + k[3][idxs] * b3Θdiff3 +
            k[4][idxs] * b4Θdiff3 + k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3 +
            k[7][idxs] * b7Θdiff3
    ) * invdt2
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing
    )
    @tsit5pre3
    @inbounds @.. broadcast = false out = (
        k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
            k[4] * b4Θdiff3 + k[5] * b5Θdiff3 + k[6] * b6Θdiff3 +
            k[7] * b7Θdiff3
    ) * invdt2
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff3 + k[2][i]*b2Θdiff3 + k[3][i]*b3Θdiff3 + k[4][i]*b4Θdiff3 + k[5][i]*b5Θdiff3 + k[6][i]*b6Θdiff3 + k[7][i]*b7Θdiff3)*invdt2
    #end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{3}}, differential_vars::Nothing
    )
    @tsit5pre3
    @views @.. broadcast = false out = (
        k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 +
            k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3 +
            k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3 +
            k[7][idxs] * b7Θdiff3
    ) * invdt2
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff3 + k[2][i]*b2Θdiff3 + k[3][i]*b3Θdiff3 + k[4][i]*b4Θdiff3 + k[5][i]*b5Θdiff3 + k[6][i]*b6Θdiff3 + k[7][i]*b7Θdiff3)*invdt2
    #end
    out
end

@def tsit5pre4 begin
    @tsit5unpack
    b1Θdiff4 = 24 * r14
    b2Θdiff4 = 24 * r24
    b3Θdiff4 = 24 * r34
    b4Θdiff4 = 24 * r44
    b5Θdiff4 = 24 * r54
    b6Θdiff4 = 24 * r64
    b7Θdiff4 = 24 * r74
    invdt3 = inv(dt)^3
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing
    )
    @tsit5pre4
    # return @.. broadcast=false k[1]*b1Θdiff4 + k[2]*b2Θdiff4 + k[3]*b3Θdiff4 + k[4]*b4Θdiff4 + k[5]*b5Θdiff4 + k[6]*b6Θdiff4 + k[7]*b7Θdiff4
    return @inbounds (
        k[1] * b1Θdiff4 + k[2] * b2Θdiff4 + k[3] * b3Θdiff4 +
            k[4] * b4Θdiff4 +
            k[5] * b5Θdiff4 + k[6] * b6Θdiff4 + k[7] * b7Θdiff4
    ) * invdt3
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{4}}, differential_vars::Nothing
    )
    @tsit5pre4
    # return @.. broadcast=false k[1][idxs]*b1Θdiff4 + k[2][idxs]*b2Θdiff4 + k[3][idxs]*b3Θdiff4 + k[4][idxs]*b4Θdiff4 + k[5][idxs]*b5Θdiff4 + k[6][idxs]*b6Θdiff4 + k[7][idxs]*b7Θdiff4
    return (
        k[1][idxs] * b1Θdiff4 + k[2][idxs] * b2Θdiff4 + k[3][idxs] * b3Θdiff4 +
            k[4][idxs] * b4Θdiff4 + k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4 +
            k[7][idxs] * b7Θdiff4
    ) * invdt3
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing
    )
    @tsit5pre4
    @inbounds @.. broadcast = false out = (
        k[1] * b1Θdiff4 + k[2] * b2Θdiff4 + k[3] * b3Θdiff4 +
            k[4] * b4Θdiff4 + k[5] * b5Θdiff4 + k[6] * b6Θdiff4 +
            k[7] * b7Θdiff4
    ) * invdt3
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff4 + k[2][i]*b2Θdiff4 + k[3][i]*b3Θdiff4 + k[4][i]*b4Θdiff4 + k[5][i]*b5Θdiff4 + k[6][i]*b6Θdiff4 + k[7][i]*b7Θdiff4)*invdt3
    #end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Tsit5ConstantCache, Tsit5CacheType}, idxs,
        T::Type{Val{4}}, differential_vars::Nothing
    )
    @tsit5pre4
    @views @.. broadcast = false out = (
        k[1][idxs] * b1Θdiff4 + k[2][idxs] * b2Θdiff4 +
            k[3][idxs] * b3Θdiff4 + k[4][idxs] * b4Θdiff4 +
            k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4 +
            k[7][idxs] * b7Θdiff4
    ) * invdt3
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff4 + k[2][i]*b2Θdiff4 + k[3][i]*b3Θdiff4 + k[4][i]*b4Θdiff4 + k[5][i]*b5Θdiff4 + k[6][i]*b6Θdiff4 + k[7][i]*b7Θdiff4)*invdt3
    #end
    out
end
