RK_WITH_SPECIAL_INTERPOLATIONS = Union{
    DP5ConstantCache, DP5Cache,
    OwrenZen3ConstantCache, OwrenZen3Cache,
    OwrenZen4ConstantCache, OwrenZen4Cache,
    OwrenZen5ConstantCache, OwrenZen5Cache,
    BS5ConstantCache, BS5Cache
}

function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::RK_WITH_SPECIAL_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::RK_WITH_SPECIAL_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars) where {D}
    throw(DerivativeOrderNotPossibleError())
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 192
"""

@def dp5pre0 begin
    b10 = Θ
    b20 = Θ * (1 - Θ)
    b30 = Θ * b20
    b40 = b20^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::DP5ConstantCache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    @dp5pre0
    @inbounds y₀ + dt * (k[1] * b10 + k[2] * b20 + k[3] * b30 + k[4] * b40)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::DP5Cache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    @dp5pre0
    @inbounds @.. broadcast=false y₀+dt *
                                     (k[1] * b10 + k[2] * b20 + k[3] * b30 + k[4] * b40)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing)
    @dp5pre0
    @views @.. broadcast=false y₀[idxs]+dt * (k[1][idxs] * b10 + k[2][idxs] * b20 +
                                         k[3][idxs] * b30 + k[4][idxs] * b40)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    @dp5pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b10 + k[2] * b20 + k[3] * b30 + k[4] * b40)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing)
    @dp5pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b10 + k[2][idxs] * b20 + k[3][idxs] * b30 +
                                    k[4][idxs] * b40)
    out
end

@def dp5pre1 begin
    b20diff = @evalpoly(Θ, 1, -2)
    b30diff = Θ * @evalpoly(Θ, 2, -3)
    b40diff = Θ * @evalpoly(Θ, 2, -6, 4)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing)
    @dp5pre1
    @inbounds @.. broadcast=false k[1]+k[2]*b20diff+k[3]*b30diff+k[4]*b40diff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    @dp5pre1
    @views @.. broadcast=false k[1][idxs]+k[2][idxs]*b20diff+k[3][idxs]*b30diff+
                               k[4][idxs]*b40diff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing)
    @dp5pre1
    @inbounds @.. broadcast=false out=k[1] + k[2] * b20diff + k[3] * b30diff +
                                      k[4] * b40diff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    @dp5pre1
    @views @.. broadcast=false out=k[1][idxs] + k[2][idxs] * b20diff +
                                   k[3][idxs] * b30diff + k[4][idxs] * b40diff
    out
end

@def dp5pre2 begin
    b20diff2 = -2
    b30diff2 = @evalpoly(Θ, 2, -6)
    b40diff2 = @evalpoly(Θ, 2, -12, 12)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{2}}, differential_vars::Nothing)
    @dp5pre2
    @inbounds @.. broadcast=false (k[2] * b20diff2 + k[3] * b30diff2 +
                                   k[4] * b40diff2)*invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing)
    @dp5pre2
    @views @.. broadcast=false (k[2][idxs] * b20diff2 + k[3][idxs] * b30diff2 +
                                k[4][idxs] * b40diff2)*invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{2}}, differential_vars::Nothing)
    @dp5pre2
    @inbounds @.. broadcast=false out=(k[2] * b20diff2 + k[3] * b30diff2 +
                                       k[4] * b40diff2) *
                                      invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{2}}, differential_vars::Nothing)
    @dp5pre2
    @views @.. broadcast=false out=(k[2][idxs] * b20diff2 + k[3][idxs] * b30diff2 +
                                    k[4][idxs] * b40diff2) * invdt
    out
end

@def dp5pre3 begin
    b30diff3 = -6
    b40diff3 = @evalpoly(Θ, -12, 24)
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{3}}, differential_vars::Nothing)
    @dp5pre3
    @inbounds @.. broadcast=false (k[3] * b30diff3 + k[4] * b40diff3)*invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{3}}, differential_vars::Nothing)
    @dp5pre3
    @views @.. broadcast=false (k[3][idxs] * b30diff3 + k[4][idxs] * b40diff3)*invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{3}}, differential_vars::Nothing)
    @dp5pre3
    @inbounds @.. broadcast=false out=(k[3] * b30diff3 + k[4] * b40diff3) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{3}}, differential_vars::Nothing)
    @dp5pre3
    @views @.. broadcast=false out=(k[3][idxs] * b30diff3 + k[4][idxs] * b40diff3) * invdt2
    out
end

@def dp5pre4 begin
    b40diff4invdt3 = 24 * inv(dt)^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{4}}, differential_vars::Nothing)
    @dp5pre4
    @inbounds @.. broadcast=false k[4]*b40diff4invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{4}}, differential_vars::Nothing)
    @dp5pre4
    @views @.. broadcast=false k[4][idxs]*b40diff4invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
        T::Type{Val{4}}, differential_vars::Nothing)
    @dp5pre4
    @inbounds @.. broadcast=false out=k[4] * b40diff4invdt3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DP5ConstantCache, DP5Cache}, idxs,
        T::Type{Val{4}}, differential_vars::Nothing)
    @dp5pre4
    @views @.. broadcast=false out=k[4][idxs] * b40diff4invdt3
    out
end

"""
"""
@def owrenzen3unpack begin
    if cache isa OrdinaryDiffEqMutableCache
        (; r13, r12, r23, r22, r33, r32) = cache.tab
    else
        (; r13, r12, r23, r22, r33, r32) = cache
    end
end

@def owrenzen3pre0 begin
    @owrenzen3unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, 1, r12, r13)
    b2Θ = Θ² * @evalpoly(Θ, r22, r23)
    b3Θ = Θ² * @evalpoly(Θ, r32, r33)
    b4Θ = Θ² * @evalpoly(Θ, -1, 1)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen3pre0
    @inbounds @.. broadcast=false y₀+dt *
                                     (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen3pre0
    @views @.. broadcast=false y₀[idxs]+dt * (k[1][idxs] * b1Θ + k[2][idxs] * b2Θ +
                                         k[3][idxs] * b3Θ +
                                         k[4][idxs] * b4Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen3pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen3pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[2][idxs] * b2Θ + k[3][idxs] * b3Θ +
                                    k[4][idxs] * b4Θ)
    out
end

@def owrenzen3pre1 begin
    @owrenzen3unpack
    b1Θdiff = @evalpoly(Θ, 1, 2*r12, 3*r13)
    b2Θdiff = Θ * @evalpoly(Θ, 2*r22, 3*r23)
    b3Θdiff = Θ * @evalpoly(Θ, 2*r32, 3*r33)
    b4Θdiff = Θ * @evalpoly(Θ, -2, 3)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen3pre1
    @inbounds @.. broadcast=false k[1]*b1Θdiff+k[2]*b2Θdiff+k[3]*b3Θdiff+k[4]*b4Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen3pre1
    @views @.. broadcast=false k[1][idxs]*b1Θdiff+k[2][idxs]*b2Θdiff+k[3][idxs]*b3Θdiff+
                               k[4][idxs]*b4Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen3pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff +
                                      k[4] * b4Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen3pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[2][idxs] * b2Θdiff +
                                   k[3][idxs] * b3Θdiff + k[4][idxs] * b4Θdiff
    out
end

@def owrenzen3pre2 begin
    @owrenzen3unpack
    b1Θdiff2 = @evalpoly(Θ, 2*r12, 6*r13)
    b2Θdiff2 = @evalpoly(Θ, 2*r22, 6*r23)
    b3Θdiff2 = @evalpoly(Θ, 2*r32, 6*r33)
    b4Θdiff2 = @evalpoly(Θ, -2, 6)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen3pre2
    @inbounds @.. broadcast=false (k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
                                   k[4] * b4Θdiff2)*invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen3pre2
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 +
                                k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2)*invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen3pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
                                       k[4] * b4Θdiff2) * invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen3pre2
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 +
                                    k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2) * invdt
    out
end

@def owrenzen3pre3 begin
    @owrenzen3unpack
    b1Θdiff3 = 6 * r13
    b2Θdiff3 = 6 * r23
    b3Θdiff3 = 6 * r33
    b4Θdiff3 = 6
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen3pre3
    @inbounds @.. broadcast=false (k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
                                   k[4] * b4Θdiff3)*invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen3pre3
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 +
                                k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3)*invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen3pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
                                       k[4] * b4Θdiff3) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
        idxs, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen3pre3
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 +
                                    k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3) * invdt2
    out
end

"""
"""
@def owrenzen4unpack begin
    if cache isa OrdinaryDiffEqMutableCache
        (; r14, r13, r12, r34, r33, r32, r44, r43, r42, r54, r53, r52, r64, r63, r62) = cache.tab
    else
        (; r14, r13, r12, r34, r33, r32, r44, r43, r42, r54, r53, r52, r64, r63, r62) = cache
    end
end

@def owrenzen4pre0 begin
    @owrenzen4unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, 1, r12, r13, r14)
    b3Θ = Θ² * @evalpoly(Θ, r32, r33, r34)
    b4Θ = Θ² * @evalpoly(Θ, r42, r43, r44)
    b5Θ = Θ² * @evalpoly(Θ, r52, r53, r54)
    b6Θ = Θ² * @evalpoly(Θ, r62, r63, r64)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen4pre0
    # return @.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen4pre0
    # return @.. broadcast=false y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[3][idxs]*b3Θ +
    #                          k[4][idxs]*b4Θ + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ)
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen4pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ +
                                       k[6] * b6Θ)
    #@inbounds for i in eachindex(out)
    #  out[i] = y₀[i] + dt*(k[1][i]*b1Θ  + k[3][i]*b3Θ + k[4][i]*b4Θ +
    #                       k[5][i]*b5Θ + k[6][i]*b6Θ)
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen4pre0
    @inbounds @.. broadcast=false out=y₀[idxs] +
                                      dt * (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ +
                                       k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
                                       k[6][idxs] * b6Θ)
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*(k[1][i]*b1Θ  + k[3][i]*b3Θ + k[4][i]*b4Θ +
    #                       k[5][i]*b5Θ + k[6][i]*b6Θ)
    #end
    out
end

@def owrenzen4pre1 begin
    @owrenzen4unpack
    b1Θdiff = @evalpoly(Θ, 1, 2*r12, 3*r13, 4*r14)
    b3Θdiff = Θ * @evalpoly(Θ, 2*r32, 3*r33, 4*r34)
    b4Θdiff = Θ * @evalpoly(Θ, 2*r42, 3*r43, 4*r44)
    b5Θdiff = Θ * @evalpoly(Θ, 2*r52, 3*r53, 4*r54)
    b6Θdiff = Θ * @evalpoly(Θ, 2*r62, 3*r63, 4*r64)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen4pre1
    @inbounds @.. broadcast=false k[1]*b1Θdiff+k[3]*b3Θdiff+k[4]*b4Θdiff+k[5]*b5Θdiff+
                                  k[6]*b6Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen4pre1
    @views @.. broadcast=false k[1][idxs]*b1Θdiff+k[3][idxs]*b3Θdiff+k[4][idxs]*b4Θdiff+
                               k[5][idxs]*b5Θdiff+k[6][idxs]*b6Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen4pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
                                      k[5] * b5Θdiff + k[6] * b6Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen4pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[3][idxs] * b3Θdiff +
                                   k[4][idxs] * b4Θdiff +
                                   k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff
    out
end

@def owrenzen4pre2 begin
    @owrenzen4unpack
    b1Θdiff2 = @evalpoly(Θ, 2*r12, 6*r13, 12*r14)
    b3Θdiff2 = @evalpoly(Θ, 2*r32, 6*r33, 12*r34)
    b4Θdiff2 = @evalpoly(Θ, 2*r42, 6*r43, 12*r44)
    b5Θdiff2 = @evalpoly(Θ, 2*r52, 6*r53, 12*r54)
    b6Θdiff2 = @evalpoly(Θ, 2*r62, 6*r63, 12*r64)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen4pre2
    @.. broadcast=false (k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                         k[5] * b5Θdiff2 + k[6] * b6Θdiff2)*invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen4pre2
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff2 + k[3][idxs] * b3Θdiff2 +
                                k[4][idxs] * b4Θdiff2 +
                                k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2)*invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen4pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                                       k[5] * b5Θdiff2 + k[6] * b6Θdiff2) * invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen4pre2
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff2 + k[3][idxs] * b3Θdiff2 +
                                    k[4][idxs] * b4Θdiff2 +
                                    k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2) * invdt
    out
end

@def owrenzen4pre3 begin
    @owrenzen4unpack
    b1Θdiff3 = @evalpoly(Θ, 6*r13, 24*r14)
    b3Θdiff3 = @evalpoly(Θ, 6*r33, 24*r34)
    b4Θdiff3 = @evalpoly(Θ, 6*r43, 24*r44)
    b5Θdiff3 = @evalpoly(Θ, 6*r53, 24*r54)
    b6Θdiff3 = @evalpoly(Θ, 6*r63, 24*r64)
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen4pre3
    @inbounds @.. broadcast=false (k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                                   k[5] * b5Θdiff3 + k[6] * b6Θdiff3)*invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen4pre3
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff3 + k[3][idxs] * b3Θdiff3 +
                                k[4][idxs] * b4Θdiff3 +
                                k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3)*invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen4pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                                       k[5] * b5Θdiff3 + k[6] * b6Θdiff3) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen4pre3
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff3 + k[3][idxs] * b3Θdiff3 +
                                    k[4][idxs] * b4Θdiff3 +
                                    k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3) * invdt2
    out
end

@def owrenzen4pre4 begin
    @owrenzen4unpack
    b1Θdiff4 = 24 * r14
    b3Θdiff4 = 24 * r34
    b4Θdiff4 = 24 * r44
    b5Θdiff4 = 24 * r54
    b6Θdiff4 = 24 * r64
    invdt3 = inv(dt)^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen4pre4
    @.. broadcast=false (k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                         k[5] * b5Θdiff4 + k[6] * b6Θdiff4)*invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen4pre4
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 +
                                k[4][idxs] * b4Θdiff4 +
                                k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4)*invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen4pre4
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                                       k[5] * b5Θdiff4 + k[6] * b6Θdiff4) * invdt3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
        idxs, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen4pre4
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 +
                                    k[4][idxs] * b4Θdiff4 +
                                    k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4) * invdt3
    out
end

"""
"""
@def owrenzen5unpack begin
    if cache isa OrdinaryDiffEqMutableCache
        (; r15, r14, r13, r12, r35, r34, r33, r32, r45, r44, r43, r42, r55, r54, r53, r52, r65, r64, r63, r62, r75, r74, r73, r72, r85, r84, r83, r82) = cache.tab
    else
        (; r15, r14, r13, r12, r35, r34, r33, r32, r45, r44, r43, r42, r55, r54, r53, r52, r65, r64, r63, r62, r75, r74, r73, r72, r85, r84, r83, r82) = cache
    end
end

@def owrenzen5pre0 begin
    @owrenzen5unpack
    Θ² = Θ * Θ
    b1Θ = Θ * @evalpoly(Θ, 1, r12, r13, r14, r15)
    b3Θ = Θ² * @evalpoly(Θ, r32, r33, r34, r35)
    b4Θ = Θ² * @evalpoly(Θ, r42, r43, r44, r45)
    b5Θ = Θ² * @evalpoly(Θ, r52, r53, r54, r55)
    b6Θ = Θ² * @evalpoly(Θ, r62, r63, r64, r65)
    b7Θ = Θ² * @evalpoly(Θ, r72, r73, r74, r75)
    b8Θ = Θ² * @evalpoly(Θ, r82, r83, r84, r85)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen5pre0
    # return @.. broadcast=false y₀ + dt*(k[1]*b1Θ  + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ +
    #                    k[7]*b7Θ + k[8]*b8Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ +
                      k[7] * b7Θ + k[8] * b8Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen5pre0
    # return @.. broadcast=false y₀[idxs] + dt*(k[1][idxs]*b1Θ  + k[3][idxs]*b3Θ +
    #                          k[4][idxs]*b4Θ + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ +
    #                          k[7][idxs]*b7Θ + k[8][idxs]*b8Θ)
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ +
            k[7][idxs] * b7Θ + k[8][idxs] * b8Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen5pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ +
                                       k[6] * b6Θ +
                                       k[7] * b7Θ + k[8] * b8Θ)
    #@inbounds for i in eachindex(out)
    #  out[i] = y₀[i] + dt*(k[1][i]*b1Θ  + k[3][i]*b3Θ + k[4][i]*b4Θ +
    #                       k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ)
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{0}}, differential_vars::Nothing)
    @owrenzen5pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt * (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ +
                                    k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ +
                                    k[7][idxs] * b7Θ + k[8][idxs] * b8Θ)
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[3][i]*b3Θ + k[4][i]*b4Θ +
    #                       k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ)
    #end
    out
end

@def owrenzen5pre1 begin
    @owrenzen5unpack
    b1Θdiff = @evalpoly(Θ, 1, 2*r12, 3*r13, 4*r14, 5*r15)
    b3Θdiff = Θ * @evalpoly(Θ, 2*r32, 3*r33, 4*r34, 5*r35)
    b4Θdiff = Θ * @evalpoly(Θ, 2*r42, 3*r43, 4*r44, 5*r45)
    b5Θdiff = Θ * @evalpoly(Θ, 2*r52, 3*r53, 4*r54, 5*r55)
    b6Θdiff = Θ * @evalpoly(Θ, 2*r62, 3*r63, 4*r64, 5*r65)
    b7Θdiff = Θ * @evalpoly(Θ, 2*r72, 3*r73, 4*r74, 5*r75)
    b8Θdiff = Θ * @evalpoly(Θ, 2*r82, 3*r83, 4*r84, 5*r85)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen5pre1
    return @inbounds k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff +
                     k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen5pre1
    k[1][idxs] * b1Θdiff + k[3][idxs] * b3Θdiff + k[4][idxs] * b4Θdiff +
    k[5][idxs] * b5Θdiff +
    k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen5pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
                                      k[5] * b5Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff +
                                      k[8] * b8Θdiff
    #@inbounds for i in eachindex(out)
    #  out[i] = k[1][i]*b1Θdiff + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff +
    #    k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff + k[8][i]*b8Θdiff
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{1}}, differential_vars::Nothing)
    @owrenzen5pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[3][idxs] * b3Θdiff +
                                   k[4][idxs] * b4Θdiff +
                                   k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
                                   k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = k[1][i]*b1Θdiff + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff +
    #    k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff + k[8][i]*b8Θdiff
    #end
    out
end

@def owrenzen5pre2 begin
    @owrenzen5unpack
    b1Θdiff2 = @evalpoly(Θ, 2*r12, 6*r13, 12*r14, 20*r15)
    b3Θdiff2 = @evalpoly(Θ, 2*r32, 6*r33, 12*r34, 20*r35)
    b4Θdiff2 = @evalpoly(Θ, 2*r42, 6*r43, 12*r44, 20*r45)
    b5Θdiff2 = @evalpoly(Θ, 2*r52, 6*r53, 12*r54, 20*r55)
    b6Θdiff2 = @evalpoly(Θ, 2*r62, 6*r63, 12*r64, 20*r65)
    b7Θdiff2 = @evalpoly(Θ, 2*r72, 6*r73, 12*r74, 20*r75)
    b8Θdiff2 = @evalpoly(Θ, 2*r82, 6*r83, 12*r84, 20*r85)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen5pre2
    return @inbounds (k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                      k[5] * b5Θdiff2 +
                      k[6] * b6Θdiff2 + k[7] * b7Θdiff2 + k[8] * b8Θdiff2) * invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen5pre2
    (k[1][idxs] * b1Θdiff2 + k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2 +
     k[5][idxs] * b5Θdiff2 +
     k[6][idxs] * b6Θdiff2 + k[7][idxs] * b7Θdiff2 + k[8][idxs] * b8Θdiff2) * invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen5pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                                       k[5] * b5Θdiff2 + k[6] * b6Θdiff2 + k[7] * b7Θdiff2 +
                                       k[8] * b8Θdiff2) * invdt
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff2 + k[3][i]*b3Θdiff2 + k[4][i]*b4Θdiff2 +
    #            k[5][i]*b5Θdiff2 + k[6][i]*b6Θdiff2 + k[7][i]*b7Θdiff2 + k[8][i]*b8Θdiff2)*invdt
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{2}}, differential_vars::Nothing)
    @owrenzen5pre2
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff2 + k[3][idxs] * b3Θdiff2 +
                                    k[4][idxs] * b4Θdiff2 +
                                    k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2 +
                                    k[7][idxs] * b7Θdiff2 + k[8][idxs] * b8Θdiff2) * invdt
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff2 + k[3][i]*b3Θdiff2 + k[4][i]*b4Θdiff2 +
    #            k[5][i]*b5Θdiff2 + k[6][i]*b6Θdiff2 + k[7][i]*b7Θdiff2 + k[8][i]*b8Θdiff2)*invdt
    #end
    out
end

@def owrenzen5pre3 begin
    @owrenzen5unpack
    b1Θdiff3 = @evalpoly(Θ, 6*r13, 24*r14, 60*r15)
    b3Θdiff3 = @evalpoly(Θ, 6*r33, 24*r34, 60*r35)
    b4Θdiff3 = @evalpoly(Θ, 6*r43, 24*r44, 60*r45)
    b5Θdiff3 = @evalpoly(Θ, 6*r53, 24*r54, 60*r55)
    b6Θdiff3 = @evalpoly(Θ, 6*r63, 24*r64, 60*r65)
    b7Θdiff3 = @evalpoly(Θ, 6*r73, 24*r74, 60*r75)
    b8Θdiff3 = @evalpoly(Θ, 6*r83, 24*r84, 60*r85)
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen5pre3
    return @inbounds (k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                      k[5] * b5Θdiff3 +
                      k[6] * b6Θdiff3 + k[7] * b7Θdiff3 + k[8] * b8Θdiff3) * invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen5pre3
    (k[1][idxs] * b1Θdiff3 + k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3 +
     k[5][idxs] * b5Θdiff3 +
     k[6][idxs] * b6Θdiff3 + k[7][idxs] * b7Θdiff3 + k[8][idxs] * b8Θdiff3) * invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen5pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                                       k[5] * b5Θdiff3 + k[6] * b6Θdiff3 + k[7] * b7Θdiff3 +
                                       k[8] * b8Θdiff3) * invdt2
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff3 + k[3][i]*b3Θdiff3 + k[4][i]*b4Θdiff3 +
    #            k[5][i]*b5Θdiff3 + k[6][i]*b6Θdiff3 + k[7][i]*b7Θdiff3 + k[8][i]*b8Θdiff3)*invdt2
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{3}}, differential_vars::Nothing)
    @owrenzen5pre3
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff3 + k[3][idxs] * b3Θdiff3 +
                                    k[4][idxs] * b4Θdiff3 +
                                    k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3 +
                                    k[7][idxs] * b7Θdiff3 + k[8][idxs] * b8Θdiff3) * invdt2
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff3 + k[3][i]*b3Θdiff3 + k[4][i]*b4Θdiff3 +
    #            k[5][i]*b5Θdiff3 + k[6][i]*b6Θdiff3 + k[7][i]*b7Θdiff3 + k[8][i]*b8Θdiff3)*invdt2
    #end
    out
end

@def owrenzen5pre4 begin
    @owrenzen5unpack
    b1Θdiff4 = @evalpoly(Θ, 24*r14, 120*r15)
    b3Θdiff4 = @evalpoly(Θ, 24*r34, 120*r35)
    b4Θdiff4 = @evalpoly(Θ, 24*r44, 120*r45)
    b5Θdiff4 = @evalpoly(Θ, 24*r54, 120*r55)
    b6Θdiff4 = @evalpoly(Θ, 24*r64, 120*r65)
    b7Θdiff4 = @evalpoly(Θ, 24*r74, 120*r75)
    b8Θdiff4 = @evalpoly(Θ, 24*r84, 120*r85)
    invdt3 = inv(dt)^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen5pre4
    return @inbounds (k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                      k[5] * b5Θdiff4 +
                      k[6] * b6Θdiff4 + k[7] * b7Θdiff4 + k[8] * b8Θdiff4) * invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen5pre4
    (k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 + k[4][idxs] * b4Θdiff4 +
     k[5][idxs] * b5Θdiff4 +
     k[6][idxs] * b6Θdiff4 + k[7][idxs] * b7Θdiff4 + k[8][idxs] * b8Θdiff4) * invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen5pre4
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                                       k[5] * b5Θdiff4 + k[6] * b6Θdiff4 + k[7] * b7Θdiff4 +
                                       k[8] * b8Θdiff4) * invdt3
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff4 + k[3][i]*b3Θdiff4 + k[4][i]*b4Θdiff4 +
    #            k[5][i]*b5Θdiff4 + k[6][i]*b6Θdiff4 + k[7][i]*b7Θdiff4 + k[8][i]*b8Θdiff4)*invdt3
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{4}}, differential_vars::Nothing)
    @owrenzen5pre4
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 +
                                    k[4][idxs] * b4Θdiff4 +
                                    k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4 +
                                    k[7][idxs] * b7Θdiff4 + k[8][idxs] * b8Θdiff4) * invdt3
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff4 + k[3][i]*b3Θdiff4 + k[4][i]*b4Θdiff4 +
    #            k[5][i]*b5Θdiff4 + k[6][i]*b6Θdiff4 + k[7][i]*b7Θdiff4 + k[8][i]*b8Θdiff4)*invdt3
    #end
    out
end

@def owrenzen5pre5 begin
    @owrenzen5unpack
    b1Θdiff5 = 120 * r15
    b3Θdiff5 = 120 * r35
    b4Θdiff5 = 120 * r45
    b5Θdiff5 = 120 * r55
    b6Θdiff5 = 120 * r65
    b7Θdiff5 = 120 * r75
    b8Θdiff5 = 120 * r85
    invdt4 = inv(dt)^4
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{5}}, differential_vars::Nothing)
    @owrenzen5pre5
    return @inbounds (k[1] * b1Θdiff5 + k[3] * b3Θdiff5 + k[4] * b4Θdiff5 +
                      k[5] * b5Θdiff5 +
                      k[6] * b6Θdiff5 + k[7] * b7Θdiff5 + k[8] * b8Θdiff5) * invdt4
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{5}}, differential_vars::Nothing)
    @owrenzen5pre5
    (k[1][idxs] * b1Θdiff5 + k[3][idxs] * b3Θdiff5 + k[4][idxs] * b4Θdiff5 +
     k[5][idxs] * b5Θdiff5 +
     k[6][idxs] * b6Θdiff5 + k[7][idxs] * b7Θdiff5 + k[8][idxs] * b8Θdiff5) * invdt4
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs::Nothing, T::Type{Val{5}}, differential_vars::Nothing)
    @owrenzen5pre5
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff5 + k[3] * b3Θdiff5 + k[4] * b4Θdiff5 +
                                       k[5] * b5Θdiff5 + k[6] * b6Θdiff5 + k[7] * b7Θdiff5 +
                                       k[8] * b8Θdiff5) * invdt4
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff5 + k[3][i]*b3Θdiff5 + k[4][i]*b4Θdiff5 +
    #            k[5][i]*b5Θdiff5 + k[6][i]*b6Θdiff5 + k[7][i]*b7Θdiff5 + k[8][i]*b8Θdiff5)*invdt4
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
        idxs, T::Type{Val{5}}, differential_vars::Nothing)
    @owrenzen5pre5
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff5 + k[3][idxs] * b3Θdiff5 +
                                    k[4][idxs] * b4Θdiff5 +
                                    k[5][idxs] * b5Θdiff5 + k[6][idxs] * b6Θdiff5 +
                                    k[7][idxs] * b7Θdiff5 + k[8][idxs] * b8Θdiff5) * invdt4
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff5 + k[3][i]*b3Θdiff5 + k[4][i]*b4Θdiff5 +
    #            k[5][i]*b5Θdiff5 + k[6][i]*b6Θdiff5 + k[7][i]*b7Θdiff5 + k[8][i]*b8Θdiff5)*invdt4
    #end
    out
end

"""
Coefficients taken from RKSuite
"""
@def bs5unpack begin
    if cache isa OrdinaryDiffEqMutableCache
        (; r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115, r114, r113, r112) = cache.tab
    else
        (; r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115, r114, r113, r112) = cache
    end
end

@def bs5pre0 begin
    @bs5unpack
    Θ² = Θ * Θ
    b1Θ = Θ² * @evalpoly(Θ, r012, r013, r014, r015, r016)
    b3Θ = Θ² * @evalpoly(Θ, r032, r033, r034, r035, r036)
    b4Θ = Θ² * @evalpoly(Θ, r042, r043, r044, r045, r046)
    b5Θ = Θ² * @evalpoly(Θ, r052, r053, r054, r055, r056)
    b6Θ = Θ² * @evalpoly(Θ, r062, r063, r064, r065, r066)
    b7Θ = Θ² * @evalpoly(Θ, r072, r073, r074, r075, r076)
    b8Θ = Θ² * @evalpoly(Θ, r082, r083, r084, r085, r086)
    b9Θ = (Θ² * Θ) * @evalpoly(Θ, r093, r094, r095,
        r096)
    b10Θ = Θ² * @evalpoly(Θ, r102, r103, r104, r105, r106)
    b11Θ = Θ² * @evalpoly(Θ, r112, r113, r114, r115, r116)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::BS5ConstantCache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    @bs5pre0
    # return @.. broadcast=false y₀ + dt*Θ*k[1] + dt*(k[1]*b1Θ  + k[3]*b3Θ + k[4]*b4Θ  + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
    return @inbounds y₀ + dt * Θ * k[1] +
                     dt * (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ +
                      k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ +
                      k[11] * b11Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::BS5Cache, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    @bs5pre0
    # return @.. broadcast=false y₀ + dt*Θ*k[1] + dt*(k[1]*b1Θ  + k[3]*b3Θ + k[4]*b4Θ  + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
    return @inbounds @.. broadcast=false y₀+dt*Θ*k[1]+
                                         dt*(k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ +
                                             k[5] * b5Θ +
                                             k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ +
                                             k[9] * b9Θ + k[10] * b10Θ + k[11] * b11Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing)
    @bs5pre0
    # return @.. broadcast=false y₀[idxs] + dt*Θ*k[1][idxs] + dt*(k[1][idxs]*b1Θ  + k[3][idxs]*b3Θ +
    #                                            k[4][idxs]*b4Θ  + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ + k[7][idxs]*b7Θ +
    #                                            k[8][idxs]*b8Θ + k[9][idxs]*b9Θ + k[10][idxs]*b10Θ + k[11][idxs]*b11Θ)
    return y₀[idxs] + dt * Θ * k[1][idxs] +
           dt * (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ +
            k[8][idxs] * b8Θ + k[9][idxs] * b9Θ + k[10][idxs] * b10Θ + k[11][idxs] * b11Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars::Nothing)
    @bs5pre0
    @inbounds @.. broadcast=false out=y₀ + dt * Θ * k[1] +
                                      dt *
                                      (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ +
                                       k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ +
                                       k[10] * b10Θ + k[11] * b11Θ)
    #@inbounds for i in eachindex(out)
    #  out[i] = y₀[i] + dt*Θ*k[1][i] + dt*(k[1][i]*b1Θ  + k[3][i]*b3Θ + k[4][i]*b4Θ  + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[10][i]*b10Θ + k[11][i]*b11Θ)
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs,
        T::Type{Val{0}}, differential_vars::Nothing)
    @bs5pre0
    @views @.. broadcast=false out=y₀[idxs] + dt * Θ * k[1][idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ + k[4][idxs] * b4Θ +
                                    k[5][idxs] * b5Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ +
                                    k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
                                    k[10][idxs] * b10Θ + k[11][idxs] * b11Θ)
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*Θ*k[1][i] + dt*(k[1][i]*b1Θ  + k[3][i]*b3Θ + k[4][i]*b4Θ  + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ + k[8][i]*b8Θ + k[9][i]*b9Θ + k[10][i]*b10Θ + k[11][i]*b11Θ)
    #end
    out
end

@def bs5pre1 begin
    @bs5unpack
    Θ² = Θ * Θ
    b1Θdiff = Θ * @evalpoly(Θ, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016)
    b3Θdiff = Θ * @evalpoly(Θ, 2*r032, 3*r033, 4*r034, 5*r035, 6*r036)
    b4Θdiff = Θ * @evalpoly(Θ, 2*r042, 3*r043, 4*r044, 5*r045, 6*r046)
    b5Θdiff = Θ * @evalpoly(Θ, 2*r052, 3*r053, 4*r054, 5*r055, 6*r056)
    b6Θdiff = Θ * @evalpoly(Θ, 2*r062, 3*r063, 4*r064, 5*r065, 6*r066)
    b7Θdiff = Θ * @evalpoly(Θ, 2*r072, 3*r073, 4*r074, 5*r075, 6*r076)
    b8Θdiff = Θ * @evalpoly(Θ, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086)
    b9Θdiff = Θ² * @evalpoly(Θ, 3*r093, 4*r094, 5*r095, 6*r096)
    b10Θdiff = Θ * @evalpoly(Θ, 2*r102, 3*r103, 4*r104, 5*r105, 6*r106)
    b11Θdiff = Θ * @evalpoly(Θ, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing)
    @bs5pre1
    # return @.. broadcast=false k[1] + k[1]*b1Θdiff  + k[3]*b3Θdiff + k[4]*b4Θdiff  + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[10]*b10Θdiff + k[11]*b11Θdiff
    return @inbounds k[1] + k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
                     k[5] * b5Θdiff +
                     k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff +
                     k[10] * b10Θdiff + k[11] * b11Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    @bs5pre1
    # return @.. broadcast=false k[1][idxs] + k[1][idxs]*b1Θdiff  + k[3][idxs]*b3Θdiff +
    #     k[4][idxs]*b4Θdiff  + k[5][idxs]*b5Θdiff + k[6][idxs]*b6Θdiff +
    #     k[7][idxs]*b7Θdiff + k[8][idxs]*b8Θdiff + k[9][idxs]*b9Θdiff +
    #     k[10][idxs]*b10Θdiff + k[11][idxs]*b11Θdiff
    return @inbounds k[1][idxs] + k[1][idxs] * b1Θdiff + k[3][idxs] * b3Θdiff +
                     k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
                     k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff + k[9][idxs] * b9Θdiff +
                     k[10][idxs] * b10Θdiff + k[11][idxs] * b11Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars::Nothing)
    @bs5pre1
    @inbounds @.. broadcast=false out=k[1] + k[1] * b1Θdiff + k[3] * b3Θdiff +
                                      k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
                                      k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff +
                                      k[10] * b10Θdiff + k[11] * b11Θdiff
    #@inbounds for i in eachindex(out)
    #  out[i] = k[1][i] + k[1][i]*b1Θdiff  + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff  + k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff + k[8][i]*b8Θdiff + k[9][i]*b9Θdiff + k[10][i]*b10Θdiff + k[11][i]*b11Θdiff
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{BS5ConstantCache, BS5Cache}, idxs,
        T::Type{Val{1}}, differential_vars::Nothing)
    @bs5pre1
    @views @.. broadcast=false out=k[1][idxs] + k[1][idxs] * b1Θdiff +
                                   k[3][idxs] * b3Θdiff + k[4][idxs] * b4Θdiff +
                                   k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
                                   k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
                                   k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff +
                                   k[11][idxs] * b11Θdiff
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = k[1][i] + k[1][i]*b1Θdiff  + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff  + k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff + k[8][i]*b8Θdiff + k[9][i]*b9Θdiff + k[10][i]*b10Θdiff + k[11][i]*b11Θdiff
    #end
    out
end
