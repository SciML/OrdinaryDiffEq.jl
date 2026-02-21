RK_WITH_SPECIAL_INTERPOLATIONS = Union{
    Vern6ConstantCache, Vern6Cache,
    Vern7ConstantCache, Vern7CacheType,
    Vern8ConstantCache, Vern8Cache,
    Vern9ConstantCache, Vern9Cache,
}

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

## Vern6
@def vern6unpack begin
    (; r011, r012, r013, r014, r015, r016, r042, r043, r044, r045, r046, r052, r053, r054, r055, r056, r062, r063, r064, r065, r066, r072, r073, r074, r075, r076, r082, r083, r084, r085, r086, r092, r093, r094, r095, r096, r102, r103, r104, r105, r106, r112, r113, r114, r115, r116, r122, r123, r124, r125, r126) = cache.tab.interp
end

@def vern6pre0 begin
    @vern6unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, r011, r012, r013, r014, r015, r016)
    b4Θ = Θ² * @evalpoly(Θ, r042, r043, r044, r045, r046)
    b5Θ = Θ² * @evalpoly(Θ, r052, r053, r054, r055, r056)
    b6Θ = Θ² * @evalpoly(Θ, r062, r063, r064, r065, r066)
    b7Θ = Θ² * @evalpoly(Θ, r072, r073, r074, r075, r076)
    b8Θ = Θ² * @evalpoly(Θ, r082, r083, r084, r085, r086)
    b9Θ = Θ² * @evalpoly(Θ, r092, r093, r094, r095, r096)
    b10Θ = Θ² * @evalpoly(Θ, r102, r103, r104, r105, r106)
    b11Θ = Θ² * @evalpoly(Θ, r112, r113, r114, r115, r116)
    b12Θ = Θ² * @evalpoly(Θ, r122, r123, r124, r125, r126)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern6ConstantCache,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern6pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ)
    return @inbounds y₀ +
        dt * (
        k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ +
            k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ + k[11] * b11Θ + k[12] * b12Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern6Cache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern6pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ)
    return @inbounds @.. broadcast = false y₀ + dt * (
        k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ +
            k[6] * b6Θ + k[7] * b7Θ +
            k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ +
            k[11] * b11Θ + k[12] * b12Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern6pre0
    return y₀[idxs] +
        dt * (
        k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
            k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
            k[10][idxs] * b10Θ + k[11][idxs] * b11Θ + k[12][idxs] * b12Θ
    )
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern6pre0
    @inbounds @.. broadcast = false out = y₀ +
        dt *
        (
        k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ +
            k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ +
            k[11] * b11Θ + k[12] * b12Θ
    )
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern6pre0
    @views @.. broadcast = false out = y₀[idxs] +
        dt *
        (
        k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
            k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ +
            k[9][idxs] * b9Θ + k[10][idxs] * b10Θ +
            k[11][idxs] * b11Θ + k[12][idxs] * b12Θ
    )
    out
end

@def vern6pre1 begin
    @vern6unpack
    b1Θdiff = @evalpoly(Θ, r011, 2 * r012, 3 * r013, 4 * r014, 5 * r015, 6 * r016)
    b4Θdiff = Θ * @evalpoly(Θ, 2 * r042, 3 * r043, 4 * r044, 5 * r045, 6 * r046)
    b5Θdiff = Θ * @evalpoly(Θ, 2 * r052, 3 * r053, 4 * r054, 5 * r055, 6 * r056)
    b6Θdiff = Θ * @evalpoly(Θ, 2 * r062, 3 * r063, 4 * r064, 5 * r065, 6 * r066)
    b7Θdiff = Θ * @evalpoly(Θ, 2 * r072, 3 * r073, 4 * r074, 5 * r075, 6 * r076)
    b8Θdiff = Θ * @evalpoly(Θ, 2 * r082, 3 * r083, 4 * r084, 5 * r085, 6 * r086)
    b9Θdiff = Θ * @evalpoly(Θ, 2 * r092, 3 * r093, 4 * r094, 5 * r095, 6 * r096)
    b10Θdiff = Θ * @evalpoly(Θ, 2 * r102, 3 * r103, 4 * r104, 5 * r105, 6 * r106)
    b11Θdiff = Θ * @evalpoly(Θ, 2 * r112, 3 * r113, 4 * r114, 5 * r115, 6 * r116)
    b12Θdiff = Θ * @evalpoly(Θ, 2 * r122, 3 * r123, 4 * r124, 5 * r125, 6 * r126)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern6pre1
    #@.. broadcast=false k[1]*b1Θdiff + k[4]*b4Θdiff + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[10]*b10Θdiff + k[11]*b11Θdiff + k[12]*b12Θdiff
    return @inbounds k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
        k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff + k[10] * b10Θdiff +
        k[11] * b11Θdiff + k[12] * b12Θdiff
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern6pre1
    return k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff +
        k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
        k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff + k[11][idxs] * b11Θdiff +
        k[12][idxs] * b12Θdiff
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern6pre1
    @inbounds @.. broadcast = false out = k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff +
        k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff +
        k[9] * b9Θdiff + k[10] * b10Θdiff + k[11] * b11Θdiff +
        k[12] * b12Θdiff
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern6pre1
    @views @.. broadcast = false out = k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff +
        k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
        k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
        k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff +
        k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff
    out
end

## Vern7
@def vern7unpack begin
    var"#T#" = constvalue(recursive_unitless_bottom_eltype(y₁))
    @OnDemandTableauExtract Vern7InterpolationCoefficients var"#T#"
end

@def vern7pre0 begin
    @vern7unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, r011, r012, r013, r014, r015, r016, r017)
    b4Θ = Θ² * @evalpoly(Θ, r042, r043, r044, r045, r046, r047)
    b5Θ = Θ² * @evalpoly(Θ, r052, r053, r054, r055, r056, r057)
    b6Θ = Θ² * @evalpoly(Θ, r062, r063, r064, r065, r066, r067)
    b7Θ = Θ² * @evalpoly(Θ, r072, r073, r074, r075, r076, r077)
    b8Θ = Θ² * @evalpoly(Θ, r082, r083, r084, r085, r086, r087)
    b9Θ = Θ² * @evalpoly(Θ, r092, r093, r094, r095, r096, r097)
    b11Θ = Θ² * @evalpoly(
        Θ, r112, r113, r114, r115, r116,
        r117
    )
    b12Θ = Θ² * @evalpoly(
        Θ, r122, r123, r124, r125, r126,
        r127
    )
    b13Θ = Θ² * @evalpoly(
        Θ, r132, r133, r134, r135, r136,
        r137
    )
    b14Θ = Θ² * @evalpoly(
        Θ, r142, r143, r144, r145, r146,
        r147
    )
    b15Θ = Θ² * @evalpoly(
        Θ, r152, r153, r154, r155, r156,
        r157
    )
    b16Θ = Θ² * @evalpoly(
        Θ, r162, r163, r164, r165, r166,
        r167
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern7ConstantCache,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern7pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[11]*b11Θ + k[12]*b12Θ + k[13]*b13Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ)
    return @inbounds y₀ +
        dt * (
        k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ +
            k[8] * b8Θ + k[9] * b9Θ + k[11] * b11Θ + k[12] * b12Θ + k[13] * b13Θ +
            k[14] * b14Θ + k[15] * b15Θ + k[16] * b16Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern7CacheType, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern7pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[11]*b11Θ + k[12]*b12Θ + k[13]*b13Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ)
    return @inbounds @.. broadcast = false y₀ + dt * (
        k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ +
            k[6] * b6Θ + k[7] * b7Θ +
            k[8] * b8Θ + k[9] * b9Θ + k[11] * b11Θ +
            k[12] * b12Θ + k[13] * b13Θ +
            k[14] * b14Θ + k[15] * b15Θ + k[16] * b16Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern7pre0
    return y₀[idxs] +
        dt * (
        k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
            k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
            k[11][idxs] * b11Θ + k[12][idxs] * b12Θ + k[13][idxs] * b13Θ +
            k[14][idxs] * b14Θ + k[15][idxs] * b15Θ + k[16][idxs] * b16Θ
    )
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern7pre0
    @inbounds @.. broadcast = false out = y₀ +
        dt *
        (
        k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ +
            k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ + k[11] * b11Θ +
            k[12] * b12Θ + k[13] * b13Θ + k[14] * b14Θ +
            k[15] * b15Θ + k[16] * b16Θ
    )
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern7pre0
    @views @.. broadcast = false out = y₀[idxs] +
        dt *
        (
        k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
            k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ +
            k[9][idxs] * b9Θ + k[11][idxs] * b11Θ +
            k[12][idxs] * b12Θ + k[13][idxs] * b13Θ +
            k[14][idxs] * b14Θ + k[15][idxs] * b15Θ +
            k[16][idxs] * b16Θ
    )
    out
end

@def vern7pre1 begin
    @vern7unpack
    b1Θdiff = @evalpoly(Θ, r011, 2 * r012, 3 * r013, 4 * r014, 5 * r015, 6 * r016, 7 * r017)
    b4Θdiff = Θ * @evalpoly(Θ, 2 * r042, 3 * r043, 4 * r044, 5 * r045, 6 * r046, 7 * r047)
    b5Θdiff = Θ * @evalpoly(Θ, 2 * r052, 3 * r053, 4 * r054, 5 * r055, 6 * r056, 7 * r057)
    b6Θdiff = Θ * @evalpoly(Θ, 2 * r062, 3 * r063, 4 * r064, 5 * r065, 6 * r066, 7 * r067)
    b7Θdiff = Θ * @evalpoly(Θ, 2 * r072, 3 * r073, 4 * r074, 5 * r075, 6 * r076, 7 * r077)
    b8Θdiff = Θ * @evalpoly(Θ, 2 * r082, 3 * r083, 4 * r084, 5 * r085, 6 * r086, 7 * r087)
    b9Θdiff = Θ * @evalpoly(Θ, 2 * r092, 3 * r093, 4 * r094, 5 * r095, 6 * r096, 7 * r097)
    b11Θdiff = Θ * @evalpoly(Θ, 2 * r112, 3 * r113, 4 * r114, 5 * r115, 6 * r116, 7 * r117)
    b12Θdiff = Θ * @evalpoly(Θ, 2 * r122, 3 * r123, 4 * r124, 5 * r125, 6 * r126, 7 * r127)
    b13Θdiff = Θ * @evalpoly(Θ, 2 * r132, 3 * r133, 4 * r134, 5 * r135, 6 * r136, 7 * r137)
    b14Θdiff = Θ * @evalpoly(Θ, 2 * r142, 3 * r143, 4 * r144, 5 * r145, 6 * r146, 7 * r147)
    b15Θdiff = Θ * @evalpoly(Θ, 2 * r152, 3 * r153, 4 * r154, 5 * r155, 6 * r156, 7 * r157)
    b16Θdiff = Θ * @evalpoly(Θ, 2 * r162, 3 * r163, 4 * r164, 5 * r165, 6 * r166, 7 * r167)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern7pre1
    #@.. broadcast=false k[1]*b1Θdiff + k[4]*b4Θdiff + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[11]*b11Θdiff + k[12]*b12Θdiff + k[13]*b13Θdiff + k[14]*b14Θdiff + k[15]*b15Θdiff + k[16]*b16Θdiff
    return @inbounds k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
        k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff + k[11] * b11Θdiff +
        k[12] * b12Θdiff + k[13] * b13Θdiff + k[14] * b14Θdiff +
        k[15] * b15Θdiff +
        k[16] * b16Θdiff
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern7pre1
    return k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff +
        k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
        k[9][idxs] * b9Θdiff + k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff +
        k[13][idxs] * b13Θdiff + k[14][idxs] * b14Θdiff + k[15][idxs] * b15Θdiff +
        k[16][idxs] * b16Θdiff
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern7pre1
    @inbounds @.. broadcast = false out = k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff +
        k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff +
        k[9] * b9Θdiff + k[11] * b11Θdiff + k[12] * b12Θdiff +
        k[13] * b13Θdiff + k[14] * b14Θdiff +
        k[15] * b15Θdiff + k[16] * b16Θdiff
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern7ConstantCache, Vern7CacheType}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern7pre1
    @views @.. broadcast = false out = k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff +
        k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
        k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
        k[9][idxs] * b9Θdiff + k[11][idxs] * b11Θdiff +
        k[12][idxs] * b12Θdiff + k[13][idxs] * b13Θdiff +
        k[14][idxs] * b14Θdiff + k[15][idxs] * b15Θdiff +
        k[16][idxs] * b16Θdiff
    out
end

## Vern8
@def vern8unpack begin
    (; r011, r012, r013, r014, r015, r016, r017, r018, r062, r063, r064, r065, r066, r067, r068, r072, r073, r074, r075, r076, r077, r078, r082, r083, r084, r085, r086, r087, r088, r092, r093, r094, r095, r096, r097, r098, r102, r103, r104, r105, r106, r107, r108, r112, r113, r114, r115, r116, r117, r118, r122, r123, r124, r125, r126, r127, r128, r142, r143, r144, r145, r146, r147, r148, r152, r153, r154, r155, r156, r157, r158, r162, r163, r164, r165, r166, r167, r168, r172, r173, r174, r175, r176, r177, r178, r182, r183, r184, r185, r186, r187, r188, r192, r193, r194, r195, r196, r197, r198, r202, r203, r204, r205, r206, r207, r208, r212, r213, r214, r215, r216, r217, r218) = cache.tab.interp
end

@def vern8pre0 begin
    @vern8unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, r011, r012, r013, r014, r015, r016, r017, r018)
    b6Θ = Θ² * @evalpoly(
        Θ, r062, r063, r064, r065, r066, r067,
        r068
    )
    b7Θ = Θ² * @evalpoly(
        Θ, r072, r073, r074, r075, r076, r077,
        r078
    )
    b8Θ = Θ² * @evalpoly(
        Θ, r082, r083, r084, r085, r086, r087,
        r088
    )
    b9Θ = Θ² * @evalpoly(
        Θ, r092, r093, r094, r095, r096, r097,
        r098
    )
    b10Θ = Θ² * @evalpoly(
        Θ, r102, r103, r104, r105, r106,
        r107, r108
    )
    b11Θ = Θ² * @evalpoly(
        Θ, r112, r113, r114, r115, r116,
        r117, r118
    )
    b12Θ = Θ² * @evalpoly(
        Θ, r122, r123, r124, r125, r126,
        r127, r128
    )
    b14Θ = Θ² * @evalpoly(
        Θ, r142, r143, r144, r145, r146,
        r147, r148
    )
    b15Θ = Θ² * @evalpoly(
        Θ, r152, r153, r154, r155, r156,
        r157, r158
    )
    b16Θ = Θ² * @evalpoly(
        Θ, r162, r163, r164, r165, r166,
        r167, r168
    )
    b17Θ = Θ² * @evalpoly(
        Θ, r172, r173, r174, r175, r176,
        r177, r178
    )
    b18Θ = Θ² * @evalpoly(
        Θ, r182, r183, r184, r185, r186,
        r187, r188
    )
    b19Θ = Θ² * @evalpoly(
        Θ, r192, r193, r194, r195, r196,
        r197, r198
    )
    b20Θ = Θ² * @evalpoly(
        Θ, r202, r203, r204, r205, r206,
        r207, r208
    )
    b21Θ = Θ² * @evalpoly(
        Θ, r212, r213, r214, r215, r216,
        r217, r218
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern8ConstantCache,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern8pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ + k[17]*b17Θ + k[18]*b18Θ + k[19]*b19Θ + k[20]*b20Θ + k[21]*b21Θ)
    return @inbounds y₀ +
        dt * (
        k[1] * b1Θ + k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ +
            k[10] * b10Θ + k[11] * b11Θ + k[12] * b12Θ + k[14] * b14Θ +
            k[15] * b15Θ +
            k[16] * b16Θ + k[17] * b17Θ + k[18] * b18Θ + k[19] * b19Θ +
            k[20] * b20Θ +
            k[21] * b21Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern8Cache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern8pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ + k[17]*b17Θ + k[18]*b18Θ + k[19]*b19Θ + k[20]*b20Θ + k[21]*b21Θ)
    return @inbounds @.. broadcast = false y₀ + dt * (
        k[1] * b1Θ + k[6] * b6Θ + k[7] * b7Θ +
            k[8] * b8Θ + k[9] * b9Θ +
            k[10] * b10Θ + k[11] * b11Θ + k[12] * b12Θ +
            k[14] * b14Θ + k[15] * b15Θ +
            k[16] * b16Θ + k[17] * b17Θ + k[18] * b18Θ +
            k[19] * b19Θ + k[20] * b20Θ +
            k[21] * b21Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern8pre0
    return y₀[idxs] +
        dt * (
        k[1][idxs] * b1Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ +
            k[8][idxs] * b8Θ + k[9][idxs] * b9Θ + k[10][idxs] * b10Θ +
            k[11][idxs] * b11Θ + k[12][idxs] * b12Θ + k[14][idxs] * b14Θ +
            k[15][idxs] * b15Θ + k[16][idxs] * b16Θ + k[17][idxs] * b17Θ +
            k[18][idxs] * b18Θ + k[19][idxs] * b19Θ + k[20][idxs] * b20Θ +
            k[21][idxs] * b21Θ
    )
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern8pre0
    @inbounds @.. broadcast = false out = y₀ +
        dt *
        (
        k[1] * b1Θ + k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ +
            k[9] * b9Θ + k[10] * b10Θ + k[11] * b11Θ +
            k[12] * b12Θ + k[14] * b14Θ + k[15] * b15Θ +
            k[16] * b16Θ + k[17] * b17Θ + k[18] * b18Θ +
            k[19] * b19Θ + k[20] * b20Θ + k[21] * b21Θ
    )
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern8pre0
    @views @.. broadcast = false out = y₀[idxs] +
        dt *
        (
        k[1][idxs] * b1Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ +
            k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
            k[10][idxs] * b10Θ + k[11][idxs] * b11Θ +
            k[12][idxs] * b12Θ + k[14][idxs] * b14Θ +
            k[15][idxs] * b15Θ + k[16][idxs] * b16Θ +
            k[17][idxs] * b17Θ + k[18][idxs] * b18Θ +
            k[19][idxs] * b19Θ + k[20][idxs] * b20Θ +
            k[21][idxs] * b21Θ
    )
    out
end

@def vern8pre1 begin
    @vern8unpack
    b1Θdiff = @evalpoly(Θ, r011, 2 * r012, 3 * r013, 4 * r014, 5 * r015, 6 * r016, 7 * r017, 8 * r018)
    b6Θdiff = Θ * @evalpoly(
        Θ, 2 * r062, 3 * r063, 4 * r064, 5 * r065, 6 * r066, 7 * r067,
        8 * r068
    )
    b7Θdiff = Θ * @evalpoly(
        Θ, 2 * r072, 3 * r073, 4 * r074, 5 * r075, 6 * r076, 7 * r077,
        8 * r078
    )
    b8Θdiff = Θ * @evalpoly(
        Θ, 2 * r082, 3 * r083, 4 * r084, 5 * r085, 6 * r086, 7 * r087,
        8 * r088
    )
    b9Θdiff = Θ * @evalpoly(
        Θ, 2 * r092, 3 * r093, 4 * r094, 5 * r095, 6 * r096, 7 * r097,
        8 * r098
    )
    b10Θdiff = Θ * @evalpoly(
        Θ, 2 * r102, 3 * r103, 4 * r104, 5 * r105, 6 * r106, 7 * r107,
        8 * r108
    )
    b11Θdiff = Θ * @evalpoly(
        Θ, 2 * r112, 3 * r113, 4 * r114, 5 * r115, 6 * r116, 7 * r117,
        8 * r118
    )
    b12Θdiff = Θ * @evalpoly(
        Θ, 2 * r122, 3 * r123, 4 * r124, 5 * r125, 6 * r126, 7 * r127,
        8 * r128
    )
    b14Θdiff = Θ * @evalpoly(
        Θ, 2 * r142, 3 * r143, 4 * r144, 5 * r145, 6 * r146, 7 * r147,
        8 * r148
    )
    b15Θdiff = Θ * @evalpoly(
        Θ, 2 * r152, 3 * r153, 4 * r154, 5 * r155, 6 * r156, 7 * r157,
        8 * r158
    )
    b16Θdiff = Θ * @evalpoly(
        Θ, 2 * r162, 3 * r163, 4 * r164, 5 * r165, 6 * r166, 7 * r167,
        8 * r168
    )
    b17Θdiff = Θ * @evalpoly(
        Θ, 2 * r172, 3 * r173, 4 * r174, 5 * r175, 6 * r176, 7 * r177,
        8 * r178
    )
    b18Θdiff = Θ * @evalpoly(
        Θ, 2 * r182, 3 * r183, 4 * r184, 5 * r185, 6 * r186, 7 * r187,
        8 * r188
    )
    b19Θdiff = Θ * @evalpoly(
        Θ, 2 * r192, 3 * r193, 4 * r194, 5 * r195, 6 * r196, 7 * r197,
        8 * r198
    )
    b20Θdiff = Θ * @evalpoly(
        Θ, 2 * r202, 3 * r203, 4 * r204, 5 * r205, 6 * r206, 7 * r207,
        8 * r208
    )
    b21Θdiff = Θ * @evalpoly(
        Θ, 2 * r212, 3 * r213, 4 * r214, 5 * r215, 6 * r216, 7 * r217,
        8 * r218
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern8pre1
    #@.. broadcast=false k[1]*b1Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[10]*b10Θdiff + k[11]*b11Θdiff + k[12]*b12Θdiff + k[14]*b14Θdiff + k[15]*b15Θdiff + k[16]*b16Θdiff + k[17]*b17Θdiff + k[18]*b18Θdiff + k[19]*b19Θdiff + k[20]*b20Θdiff + k[21]*b21Θdiff
    return @inbounds k[1] * b1Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff +
        k[9] * b9Θdiff + k[10] * b10Θdiff + k[11] * b11Θdiff +
        k[12] * b12Θdiff +
        k[14] * b14Θdiff + k[15] * b15Θdiff + k[16] * b16Θdiff +
        k[17] * b17Θdiff +
        k[18] * b18Θdiff + k[19] * b19Θdiff + k[20] * b20Θdiff +
        k[21] * b21Θdiff
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern8pre1
    return k[1][idxs] * b1Θdiff + k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff +
        k[8][idxs] * b8Θdiff + k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff +
        k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff + k[14][idxs] * b14Θdiff +
        k[15][idxs] * b15Θdiff + k[16][idxs] * b16Θdiff + k[17][idxs] * b17Θdiff +
        k[18][idxs] * b18Θdiff + k[19][idxs] * b19Θdiff + k[20][idxs] * b20Θdiff +
        k[21][idxs] * b21Θdiff
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern8pre1
    @inbounds @.. broadcast = false out = k[1] * b1Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff +
        k[8] * b8Θdiff + k[9] * b9Θdiff + k[10] * b10Θdiff +
        k[11] * b11Θdiff + k[12] * b12Θdiff +
        k[14] * b14Θdiff + k[15] * b15Θdiff +
        k[16] * b16Θdiff + k[17] * b17Θdiff +
        k[18] * b18Θdiff + k[19] * b19Θdiff +
        k[20] * b20Θdiff + k[21] * b21Θdiff
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern8pre1
    @views @.. broadcast = false out = k[1][idxs] * b1Θdiff + k[6][idxs] * b6Θdiff +
        k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
        k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff +
        k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff +
        k[14][idxs] * b14Θdiff + k[15][idxs] * b15Θdiff +
        k[16][idxs] * b16Θdiff + k[17][idxs] * b17Θdiff +
        k[18][idxs] * b18Θdiff + k[19][idxs] * b19Θdiff +
        k[20][idxs] * b20Θdiff + k[21][idxs] * b21Θdiff
    out
end

## Vern9
@def vern9unpack begin
    var"#T#" = constvalue(recursive_unitless_bottom_eltype(y₁))
    @OnDemandTableauExtract Vern9InterpolationCoefficients var"#T#"
end

@def vern9pre0 begin
    @vern9unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(
        Θ, r011, r012, r013, r014, r015, r016, r017, r018,
        r019
    )
    b8Θ = Θ² * @evalpoly(
        Θ, r082, r083, r084, r085, r086, r087,
        r088, r089
    )
    b9Θ = Θ² * @evalpoly(
        Θ, r092, r093, r094, r095, r096, r097,
        r098, r099
    )
    b10Θ = Θ² * @evalpoly(
        Θ, r102, r103, r104, r105, r106,
        r107, r108, r109
    )
    b11Θ = Θ² * @evalpoly(
        Θ, r112, r113, r114, r115, r116,
        r117, r118, r119
    )
    b12Θ = Θ² * @evalpoly(
        Θ, r122, r123, r124, r125, r126,
        r127, r128, r129
    )
    b13Θ = Θ² * @evalpoly(
        Θ, r132, r133, r134, r135, r136,
        r137, r138, r139
    )
    b14Θ = Θ² * @evalpoly(
        Θ, r142, r143, r144, r145, r146,
        r147, r148, r149
    )
    b15Θ = Θ² * @evalpoly(
        Θ, r152, r153, r154, r155, r156,
        r157, r158, r159
    )
    b17Θ = Θ² * @evalpoly(
        Θ, r172, r173, r174, r175, r176,
        r177, r178, r179
    )
    b18Θ = Θ² * @evalpoly(
        Θ, r182, r183, r184, r185, r186,
        r187, r188, r189
    )
    b19Θ = Θ² * @evalpoly(
        Θ, r192, r193, r194, r195, r196,
        r197, r198, r199
    )
    b20Θ = Θ² * @evalpoly(
        Θ, r202, r203, r204, r205, r206,
        r207, r208, r209
    )
    b21Θ = Θ² * @evalpoly(
        Θ, r212, r213, r214, r215, r216,
        r217, r218, r219
    )
    b22Θ = Θ² * @evalpoly(
        Θ, r222, r223, r224, r225, r226,
        r227, r228, r229
    )
    b23Θ = Θ² * @evalpoly(
        Θ, r232, r233, r234, r235, r236,
        r237, r238, r239
    )
    b24Θ = Θ² * @evalpoly(
        Θ, r242, r243, r244, r245, r246,
        r247, r248, r249
    )
    b25Θ = Θ² * @evalpoly(
        Θ, r252, r253, r254, r255, r256,
        r257, r258, r259
    )
    b26Θ = Θ² * @evalpoly(
        Θ, r262, r263, r264, r265, r266,
        r267, r268, r269
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern9ConstantCache,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern9pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[2]*b8Θ + k[3]*b9Θ + k[4]*b10Θ + k[5]*b11Θ + k[6]*b12Θ + k[7]*b13Θ + k[8]*b14Θ + k[9]*b15Θ + k[11]*b17Θ + k[12]*b18Θ + k[13]*b19Θ + k[14]*b20Θ + k[15]*b21Θ + k[16]*b22Θ + k[17]*b23Θ + k[18]*b24Θ + k[19]*b25Θ + k[20]*b26Θ)
    return @inbounds y₀ +
        dt *
        (
        k[1] * b1Θ + k[2] * b8Θ + k[3] * b9Θ + k[4] * b10Θ + k[5] * b11Θ +
            k[6] * b12Θ + k[7] * b13Θ + k[8] * b14Θ + k[9] * b15Θ + k[11] * b17Θ +
            k[12] * b18Θ + k[13] * b19Θ + k[14] * b20Θ + k[15] * b21Θ +
            k[16] * b22Θ +
            k[17] * b23Θ + k[18] * b24Θ + k[19] * b25Θ + k[20] * b26Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::Vern9Cache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern9pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[2]*b8Θ + k[3]*b9Θ + k[4]*b10Θ + k[5]*b11Θ + k[6]*b12Θ + k[7]*b13Θ + k[8]*b14Θ + k[9]*b15Θ + k[11]*b17Θ + k[12]*b18Θ + k[13]*b19Θ + k[14]*b20Θ + k[15]*b21Θ + k[16]*b22Θ + k[17]*b23Θ + k[18]*b24Θ + k[19]*b25Θ + k[20]*b26Θ)
    return @inbounds @.. broadcast = false y₀ + dt * (
        k[1] * b1Θ + k[2] * b8Θ + k[3] * b9Θ +
            k[4] * b10Θ + k[5] * b11Θ +
            k[6] * b12Θ + k[7] * b13Θ + k[8] * b14Θ +
            k[9] * b15Θ + k[11] * b17Θ +
            k[12] * b18Θ + k[13] * b19Θ + k[14] * b20Θ +
            k[15] * b21Θ + k[16] * b22Θ +
            k[17] * b23Θ + k[18] * b24Θ + k[19] * b25Θ +
            k[20] * b26Θ
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern9pre0
    return y₀[idxs] +
        dt * (
        k[1][idxs] * b1Θ + k[2][idxs] * b8Θ + k[3][idxs] * b9Θ +
            k[4][idxs] * b10Θ + k[5][idxs] * b11Θ + k[6][idxs] * b12Θ +
            k[7][idxs] * b13Θ + k[8][idxs] * b14Θ + k[9][idxs] * b15Θ + k[11][idxs] * b17Θ +
            k[12][idxs] * b18Θ + k[13][idxs] * b19Θ + k[14][idxs] * b20Θ +
            k[15][idxs] * b21Θ +
            k[16][idxs] * b22Θ + k[17][idxs] * b23Θ + k[18][idxs] * b24Θ +
            k[19][idxs] * b25Θ + k[20][idxs] * b26Θ
    )
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern9pre0
    @inbounds @.. broadcast = false out = y₀ +
        dt *
        (
        k[1] * b1Θ + k[2] * b8Θ + k[3] * b9Θ + k[4] * b10Θ +
            k[5] * b11Θ + k[6] * b12Θ + k[7] * b13Θ +
            k[8] * b14Θ + k[9] * b15Θ + k[11] * b17Θ +
            k[12] * b18Θ + k[13] * b19Θ + k[14] * b20Θ +
            k[15] * b21Θ + k[16] * b22Θ + k[17] * b23Θ +
            k[18] * b24Θ + k[19] * b25Θ + k[20] * b26Θ
    )
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @vern9pre0
    @views @.. broadcast = false out = y₀[idxs] +
        dt *
        (
        k[1][idxs] * b1Θ + k[2][idxs] * b8Θ + k[3][idxs] * b9Θ +
            k[4][idxs] * b10Θ + k[5][idxs] * b11Θ +
            k[6][idxs] * b12Θ + k[7][idxs] * b13Θ +
            k[8][idxs] * b14Θ + k[9][idxs] * b15Θ +
            k[11][idxs] * b17Θ + k[12][idxs] * b18Θ +
            k[13][idxs] * b19Θ + k[14][idxs] * b20Θ +
            k[15][idxs] * b21Θ + k[16][idxs] * b22Θ +
            k[17][idxs] * b23Θ + k[18][idxs] * b24Θ +
            k[19][idxs] * b25Θ + k[20][idxs] * b26Θ
    )
    out
end

@def vern9pre1 begin
    @vern9unpack
    b1Θdiff = @evalpoly(
        Θ, r011, 2 * r012, 3 * r013, 4 * r014, 5 * r015, 6 * r016, 7 * r017, 8 * r018,
        9 * r019
    )
    b8Θdiff = Θ * @evalpoly(
        Θ, 2 * r082, 3 * r083, 4 * r084, 5 * r085, 6 * r086, 7 * r087,
        8 * r088,
        9 * r089
    )
    b9Θdiff = Θ * @evalpoly(
        Θ, 2 * r092, 3 * r093, 4 * r094, 5 * r095, 6 * r096, 7 * r097,
        8 * r098,
        9 * r099
    )
    b10Θdiff = Θ * @evalpoly(
        Θ, 2 * r102, 3 * r103, 4 * r104, 5 * r105, 6 * r106, 7 * r107,
        8 * r108,
        9 * r109
    )
    b11Θdiff = Θ * @evalpoly(
        Θ, 2 * r112, 3 * r113, 4 * r114, 5 * r115, 6 * r116, 7 * r117,
        8 * r118,
        9 * r119
    )
    b12Θdiff = Θ * @evalpoly(
        Θ, 2 * r122, 3 * r123, 4 * r124, 5 * r125, 6 * r126, 7 * r127,
        8 * r128,
        9 * r129
    )
    b13Θdiff = Θ * @evalpoly(
        Θ, 2 * r132, 3 * r133, 4 * r134, 5 * r135, 6 * r136, 7 * r137,
        8 * r138,
        9 * r139
    )
    b14Θdiff = Θ * @evalpoly(
        Θ, 2 * r142, 3 * r143, 4 * r144, 5 * r145, 6 * r146, 7 * r147,
        8 * r148,
        9 * r149
    )
    b15Θdiff = Θ * @evalpoly(
        Θ, 2 * r152, 3 * r153, 4 * r154, 5 * r155, 6 * r156, 7 * r157,
        8 * r158,
        9 * r159
    )
    b17Θdiff = Θ * @evalpoly(
        Θ, 2 * r172, 3 * r173, 4 * r174, 5 * r175, 6 * r176, 7 * r177,
        8 * r178,
        9 * r179
    )
    b18Θdiff = Θ * @evalpoly(
        Θ, 2 * r182, 3 * r183, 4 * r184, 5 * r185, 6 * r186, 7 * r187,
        8 * r188,
        9 * r189
    )
    b19Θdiff = Θ * @evalpoly(
        Θ, 2 * r192, 3 * r193, 4 * r194, 5 * r195, 6 * r196, 7 * r197,
        8 * r198,
        9 * r199
    )
    b20Θdiff = Θ * @evalpoly(
        Θ, 2 * r202, 3 * r203, 4 * r204, 5 * r205, 6 * r206, 7 * r207,
        8 * r208,
        9 * r209
    )
    b21Θdiff = Θ * @evalpoly(
        Θ, 2 * r212, 3 * r213, 4 * r214, 5 * r215, 6 * r216, 7 * r217,
        8 * r218,
        9 * r219
    )
    b22Θdiff = Θ * @evalpoly(
        Θ, 2 * r222, 3 * r223, 4 * r224, 5 * r225, 6 * r226, 7 * r227,
        8 * r228,
        9 * r229
    )
    b23Θdiff = Θ * @evalpoly(
        Θ, 2 * r232, 3 * r233, 4 * r234, 5 * r235, 6 * r236, 7 * r237,
        8 * r238,
        9 * r239
    )
    b24Θdiff = Θ * @evalpoly(
        Θ, 2 * r242, 3 * r243, 4 * r244, 5 * r245, 6 * r246, 7 * r247,
        8 * r248,
        9 * r249
    )
    b25Θdiff = Θ * @evalpoly(
        Θ, 2 * r252, 3 * r253, 4 * r254, 5 * r255, 6 * r256, 7 * r257,
        8 * r258,
        9 * r259
    )
    b26Θdiff = Θ * @evalpoly(
        Θ, 2 * r262, 3 * r263, 4 * r264, 5 * r265, 6 * r266, 7 * r267,
        8 * r268,
        9 * r269
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern9pre1
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
        k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff + k[8] * b14Θdiff +
        k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
        k[13] * b19Θdiff +
        k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
        k[17] * b23Θdiff +
        k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern9pre1
    return k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
        k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
        k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
        k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
        k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
        k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
        k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
        k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
        k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern9pre1
    @inbounds @.. broadcast = false out = k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
        k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
        k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
        k[11] * b17Θdiff + k[12] * b18Θdiff +
        k[13] * b19Θdiff + k[14] * b20Θdiff +
        k[15] * b21Θdiff + k[16] * b22Θdiff +
        k[17] * b23Θdiff + k[18] * b24Θdiff +
        k[19] * b25Θdiff + k[20] * b26Θdiff
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing
    )
    @vern9pre1
    @views @.. broadcast = false out = k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
        k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
        k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
        k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
        k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
        k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
        k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
        k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
        k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
        k[20][idxs] * b26Θdiff
    out
end

@def vern9pre2 begin
    @vern9unpack
    b1Θdiff = @evalpoly(
        Θ, 2 * r012, 6 * r013, 12 * r014, 20 * r015, 30 * r016, 42 * r017, 56 * r018,
        72 * r019
    )
    b8Θdiff = @evalpoly(
        Θ, 2 * r082, 6 * r083, 12 * r084, 20 * r085, 30 * r086, 42 * r087, 56 * r088,
        72 * r089
    )
    b9Θdiff = @evalpoly(
        Θ, 2 * r092, 6 * r093, 12 * r094, 20 * r095, 30 * r096, 42 * r097, 56 * r098,
        72 * r099
    )
    b10Θdiff = @evalpoly(
        Θ, 2 * r102, 6 * r103, 12 * r104, 20 * r105, 30 * r106, 42 * r107, 56 * r108,
        72 * r109
    )
    b11Θdiff = @evalpoly(
        Θ, 2 * r112, 6 * r113, 12 * r114, 20 * r115, 30 * r116, 42 * r117, 56 * r118,
        72 * r119
    )
    b12Θdiff = @evalpoly(
        Θ, 2 * r122, 6 * r123, 12 * r124, 20 * r125, 30 * r126, 42 * r127, 56 * r128,
        72 * r129
    )
    b13Θdiff = @evalpoly(
        Θ, 2 * r132, 6 * r133, 12 * r134, 20 * r135, 30 * r136, 42 * r137, 56 * r138,
        72 * r139
    )
    b14Θdiff = @evalpoly(
        Θ, 2 * r142, 6 * r143, 12 * r144, 20 * r145, 30 * r146, 42 * r147, 56 * r148,
        72 * r149
    )
    b15Θdiff = @evalpoly(
        Θ, 2 * r152, 6 * r153, 12 * r154, 20 * r155, 30 * r156, 42 * r157, 56 * r158,
        72 * r159
    )
    b17Θdiff = @evalpoly(
        Θ, 2 * r172, 6 * r173, 12 * r174, 20 * r175, 30 * r176, 42 * r177, 56 * r178,
        72 * r179
    )
    b18Θdiff = @evalpoly(
        Θ, 2 * r182, 6 * r183, 12 * r184, 20 * r185, 30 * r186, 42 * r187, 56 * r188,
        72 * r189
    )
    b19Θdiff = @evalpoly(
        Θ, 2 * r192, 6 * r193, 12 * r194, 20 * r195, 30 * r196, 42 * r197, 56 * r198,
        72 * r199
    )
    b20Θdiff = @evalpoly(
        Θ, 2 * r202, 6 * r203, 12 * r204, 20 * r205, 30 * r206, 42 * r207, 56 * r208,
        72 * r209
    )
    b21Θdiff = @evalpoly(
        Θ, 2 * r212, 6 * r213, 12 * r214, 20 * r215, 30 * r216, 42 * r217, 56 * r218,
        72 * r219
    )
    b22Θdiff = @evalpoly(
        Θ, 2 * r222, 6 * r223, 12 * r224, 20 * r225, 30 * r226, 42 * r227, 56 * r228,
        72 * r229
    )
    b23Θdiff = @evalpoly(
        Θ, 2 * r232, 6 * r233, 12 * r234, 20 * r235, 30 * r236, 42 * r237, 56 * r238,
        72 * r239
    )
    b24Θdiff = @evalpoly(
        Θ, 2 * r242, 6 * r243, 12 * r244, 20 * r245, 30 * r246, 42 * r247, 56 * r248,
        72 * r249
    )
    b25Θdiff = @evalpoly(
        Θ, 2 * r252, 6 * r253, 12 * r254, 20 * r255, 30 * r256, 42 * r257, 56 * r258,
        72 * r259
    )
    b26Θdiff = @evalpoly(
        Θ, 2 * r262, 6 * r263, 12 * r264, 20 * r265, 30 * r266, 42 * r267, 56 * r268,
        72 * r269
    )
    invdt = inv(dt)
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing
    )
    @vern9pre2
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds (
        k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
            k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff +
            k[8] * b14Θdiff +
            k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
            k[13] * b19Θdiff +
            k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
            k[17] * b23Θdiff +
            k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff
    ) * invdt
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing
    )
    @vern9pre2
    return (
        k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
            k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
            k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
            k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
            k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
            k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
            k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
            k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
            k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff
    ) * invdt
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing
    )
    @vern9pre2
    @inbounds @.. broadcast = false out = (
        k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
            k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
            k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
            k[11] * b17Θdiff + k[12] * b18Θdiff +
            k[13] * b19Θdiff + k[14] * b20Θdiff +
            k[15] * b21Θdiff + k[16] * b22Θdiff +
            k[17] * b23Θdiff + k[18] * b24Θdiff +
            k[19] * b25Θdiff + k[20] * b26Θdiff
    ) * invdt
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing
    )
    @vern9pre2
    @views @.. broadcast = false out = (
        k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
            k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
            k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
            k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
            k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
            k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
            k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
            k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
            k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
            k[20][idxs] * b26Θdiff
    ) * invdt
    out
end

@def vern9pre3 begin
    @vern9unpack
    b1Θdiff = @evalpoly(Θ, 6 * r013, 24 * r014, 60 * r015, 120 * r016, 210 * r017, 336 * r018, 504 * r019)
    b8Θdiff = @evalpoly(Θ, 6 * r083, 24 * r084, 60 * r085, 120 * r086, 210 * r087, 336 * r088, 504 * r089)
    b9Θdiff = @evalpoly(Θ, 6 * r093, 24 * r094, 60 * r095, 120 * r096, 210 * r097, 336 * r098, 504 * r099)
    b10Θdiff = @evalpoly(
        Θ, 6 * r103, 24 * r104, 60 * r105, 120 * r106, 210 * r107, 336 * r108,
        504 * r109
    )
    b11Θdiff = @evalpoly(
        Θ, 6 * r113, 24 * r114, 60 * r115, 120 * r116, 210 * r117, 336 * r118,
        504 * r119
    )
    b12Θdiff = @evalpoly(
        Θ, 6 * r123, 24 * r124, 60 * r125, 120 * r126, 210 * r127, 336 * r128,
        504 * r129
    )
    b13Θdiff = @evalpoly(
        Θ, 6 * r133, 24 * r134, 60 * r135, 120 * r136, 210 * r137, 336 * r138,
        504 * r139
    )
    b14Θdiff = @evalpoly(
        Θ, 6 * r143, 24 * r144, 60 * r145, 120 * r146, 210 * r147, 336 * r148,
        504 * r149
    )
    b15Θdiff = @evalpoly(
        Θ, 6 * r153, 24 * r154, 60 * r155, 120 * r156, 210 * r157, 336 * r158,
        504 * r159
    )
    b17Θdiff = @evalpoly(
        Θ, 6 * r173, 24 * r174, 60 * r175, 120 * r176, 210 * r177, 336 * r178,
        504 * r179
    )
    b18Θdiff = @evalpoly(
        Θ, 6 * r183, 24 * r184, 60 * r185, 120 * r186, 210 * r187, 336 * r188,
        504 * r189
    )
    b19Θdiff = @evalpoly(
        Θ, 6 * r193, 24 * r194, 60 * r195, 120 * r196, 210 * r197, 336 * r198,
        504 * r199
    )
    b20Θdiff = @evalpoly(
        Θ, 6 * r203, 24 * r204, 60 * r205, 120 * r206, 210 * r207, 336 * r208,
        504 * r209
    )
    b21Θdiff = @evalpoly(
        Θ, 6 * r213, 24 * r214, 60 * r215, 120 * r216, 210 * r217, 336 * r218,
        504 * r219
    )
    b22Θdiff = @evalpoly(
        Θ, 6 * r223, 24 * r224, 60 * r225, 120 * r226, 210 * r227, 336 * r228,
        504 * r229
    )
    b23Θdiff = @evalpoly(
        Θ, 6 * r233, 24 * r234, 60 * r235, 120 * r236, 210 * r237, 336 * r238,
        504 * r239
    )
    b24Θdiff = @evalpoly(
        Θ, 6 * r243, 24 * r244, 60 * r245, 120 * r246, 210 * r247, 336 * r248,
        504 * r249
    )
    b25Θdiff = @evalpoly(
        Θ, 6 * r253, 24 * r254, 60 * r255, 120 * r256, 210 * r257, 336 * r258,
        504 * r259
    )
    b26Θdiff = @evalpoly(
        Θ, 6 * r263, 24 * r264, 60 * r265, 120 * r266, 210 * r267, 336 * r268,
        504 * r269
    )
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing
    )
    @vern9pre3
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds (
        k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
            k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff +
            k[8] * b14Θdiff +
            k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
            k[13] * b19Θdiff +
            k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
            k[17] * b23Θdiff +
            k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff
    ) * invdt2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{3}}, differential_vars::Nothing
    )
    @vern9pre3
    return (
        k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
            k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
            k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
            k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
            k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
            k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
            k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
            k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
            k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff
    ) * invdt2
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing
    )
    @vern9pre3
    @inbounds @.. broadcast = false out = (
        k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
            k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
            k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
            k[11] * b17Θdiff + k[12] * b18Θdiff +
            k[13] * b19Θdiff + k[14] * b20Θdiff +
            k[15] * b21Θdiff + k[16] * b22Θdiff +
            k[17] * b23Θdiff + k[18] * b24Θdiff +
            k[19] * b25Θdiff + k[20] * b26Θdiff
    ) * invdt2
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{3}}, differential_vars::Nothing
    )
    @vern9pre3
    @views @.. broadcast = false out = (
        k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
            k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
            k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
            k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
            k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
            k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
            k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
            k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
            k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
            k[20][idxs] * b26Θdiff
    ) * invdt2
    out
end

@def vern9pre4 begin
    @vern9unpack
    b1Θdiff = @evalpoly(Θ, 24 * r014, 120 * r015, 360 * r016, 840 * r017, 1680 * r018, 3024 * r019)
    b8Θdiff = @evalpoly(Θ, 24 * r084, 120 * r085, 360 * r086, 840 * r087, 1680 * r088, 3024 * r089)
    b9Θdiff = @evalpoly(Θ, 24 * r094, 120 * r095, 360 * r096, 840 * r097, 1680 * r098, 3024 * r099)
    b10Θdiff = @evalpoly(Θ, 24 * r104, 120 * r105, 360 * r106, 840 * r107, 1680 * r108, 3024 * r109)
    b11Θdiff = @evalpoly(Θ, 24 * r114, 120 * r115, 360 * r116, 840 * r117, 1680 * r118, 3024 * r119)
    b12Θdiff = @evalpoly(Θ, 24 * r124, 120 * r125, 360 * r126, 840 * r127, 1680 * r128, 3024 * r129)
    b13Θdiff = @evalpoly(Θ, 24 * r134, 120 * r135, 360 * r136, 840 * r137, 1680 * r138, 3024 * r139)
    b14Θdiff = @evalpoly(Θ, 24 * r144, 120 * r145, 360 * r146, 840 * r147, 1680 * r148, 3024 * r149)
    b15Θdiff = @evalpoly(Θ, 24 * r154, 120 * r155, 360 * r156, 840 * r157, 1680 * r158, 3024 * r159)
    b17Θdiff = @evalpoly(Θ, 24 * r174, 120 * r175, 360 * r176, 840 * r177, 1680 * r178, 3024 * r179)
    b18Θdiff = @evalpoly(Θ, 24 * r184, 120 * r185, 360 * r186, 840 * r187, 1680 * r188, 3024 * r189)
    b19Θdiff = @evalpoly(Θ, 24 * r194, 120 * r195, 360 * r196, 840 * r197, 1680 * r198, 3024 * r199)
    b20Θdiff = @evalpoly(Θ, 24 * r204, 120 * r205, 360 * r206, 840 * r207, 1680 * r208, 3024 * r209)
    b21Θdiff = @evalpoly(Θ, 24 * r214, 120 * r215, 360 * r216, 840 * r217, 1680 * r218, 3024 * r219)
    b22Θdiff = @evalpoly(Θ, 24 * r224, 120 * r225, 360 * r226, 840 * r227, 1680 * r228, 3024 * r229)
    b23Θdiff = @evalpoly(Θ, 24 * r234, 120 * r235, 360 * r236, 840 * r237, 1680 * r238, 3024 * r239)
    b24Θdiff = @evalpoly(Θ, 24 * r244, 120 * r245, 360 * r246, 840 * r247, 1680 * r248, 3024 * r249)
    b25Θdiff = @evalpoly(Θ, 24 * r254, 120 * r255, 360 * r256, 840 * r257, 1680 * r258, 3024 * r259)
    b26Θdiff = @evalpoly(Θ, 24 * r264, 120 * r265, 360 * r266, 840 * r267, 1680 * r268, 3024 * r269)
    invdt3 = inv(dt)^3
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing
    )
    @vern9pre4
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds (
        k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
            k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff +
            k[8] * b14Θdiff +
            k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
            k[13] * b19Θdiff +
            k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
            k[17] * b23Θdiff +
            k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff
    ) * invdt3
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{4}}, differential_vars::Nothing
    )
    @vern9pre4
    return (
        k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
            k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
            k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
            k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
            k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
            k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
            k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
            k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
            k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff
    ) * invdt3
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing
    )
    @vern9pre4
    @inbounds @.. broadcast = false out = (
        k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
            k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
            k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
            k[11] * b17Θdiff + k[12] * b18Θdiff +
            k[13] * b19Θdiff + k[14] * b20Θdiff +
            k[15] * b21Θdiff + k[16] * b22Θdiff +
            k[17] * b23Θdiff + k[18] * b24Θdiff +
            k[19] * b25Θdiff + k[20] * b26Θdiff
    ) * invdt3
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
        T::Type{Val{4}}, differential_vars::Nothing
    )
    @vern9pre4
    @views @.. broadcast = false out = (
        k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
            k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
            k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
            k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
            k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
            k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
            k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
            k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
            k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
            k[20][idxs] * b26Θdiff
    ) * invdt3
    out
end
