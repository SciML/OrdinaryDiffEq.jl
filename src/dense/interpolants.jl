@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{FunctionMapConstantCache, FunctionMapCache},
                                  idxs, T::Type{Val{0}})
    y₀
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{FunctionMapConstantCache, FunctionMapCache},
                                   idxs, T::Type{Val{0}})
    recursivecopy!(out, y₀)
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 192
"""
@def dp5pre0 begin
    b10 = @evalpoly(Θ, 0, 1)
    b20 = @evalpoly(Θ, 0, 1, -1)
    b30 = @evalpoly(Θ, 0, 0, 1, -1)
    b40 = @evalpoly(Θ, 0, 0, 1, -2, 1)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::DP5ConstantCache, idxs::Nothing,
                                  T::Type{Val{0}})
    @dp5pre0
    @inbounds y₀ + dt * (k[1] * b10 + k[2] * b20 + k[3] * b30 + k[4] * b40)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::DP5Cache, idxs::Nothing,
                                  T::Type{Val{0}})
    @dp5pre0
    @inbounds @.. broadcast=false y₀+dt *
                                     (k[1] * b10 + k[2] * b20 + k[3] * b30 + k[4] * b40)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                  T::Type{Val{0}})
    @dp5pre0
    @views @.. broadcast=false y₀[idxs]+dt * (k[1][idxs] * b10 + k[2][idxs] * b20 +
                                         k[3][idxs] * b30 + k[4][idxs] * b40)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                   T::Type{Val{0}})
    @dp5pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b10 + k[2] * b20 + k[3] * b30 + k[4] * b40)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                   T::Type{Val{0}})
    @dp5pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b10 + k[2][idxs] * b20 + k[3][idxs] * b30 +
                                    k[4][idxs] * b40)
    out
end

@def dp5pre1 begin
    b20diff = @evalpoly(Θ, 1, -2)
    b30diff = @evalpoly(Θ, 0, 2, -3)
    b40diff = @evalpoly(Θ, 0, 2, -6, 4)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                  T::Type{Val{1}})
    @dp5pre1
    @inbounds @.. broadcast=false k[1]+k[2]*b20diff+k[3]*b30diff+k[4]*b40diff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                  T::Type{Val{1}})
    @dp5pre1
    @views @.. broadcast=false k[1][idxs]+k[2][idxs]*b20diff+k[3][idxs]*b30diff+
                               k[4][idxs]*b40diff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                   T::Type{Val{1}})
    @dp5pre1
    @inbounds @.. broadcast=false out=k[1] + k[2] * b20diff + k[3] * b30diff +
                                      k[4] * b40diff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                   T::Type{Val{1}})
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
                                  T::Type{Val{2}})
    @dp5pre2
    @inbounds @.. broadcast=false (k[2] * b20diff2 + k[3] * b30diff2 + k[4] * b40diff2)*invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                  T::Type{Val{2}})
    @dp5pre2
    @views @.. broadcast=false (k[2][idxs] * b20diff2 + k[3][idxs] * b30diff2 +
                                k[4][idxs] * b40diff2)*invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                   T::Type{Val{2}})
    @dp5pre2
    @inbounds @.. broadcast=false out=(k[2] * b20diff2 + k[3] * b30diff2 + k[4] * b40diff2) *
                                      invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                   T::Type{Val{2}})
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
                                  T::Type{Val{3}})
    @dp5pre3
    @inbounds @.. broadcast=false (k[3] * b30diff3 + k[4] * b40diff3)*invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                  T::Type{Val{3}})
    @dp5pre3
    @views @.. broadcast=false (k[3][idxs] * b30diff3 + k[4][idxs] * b40diff3)*invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                   T::Type{Val{3}})
    @dp5pre3
    @inbounds @.. broadcast=false out=(k[3] * b30diff3 + k[4] * b40diff3) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                   T::Type{Val{3}})
    @dp5pre3
    @views @.. broadcast=false out=(k[3][idxs] * b30diff3 + k[4][idxs] * b40diff3) * invdt2
    out
end

@def dp5pre4 begin b40diff4invdt3 = 24 * inv(dt)^3 end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                  T::Type{Val{4}})
    @dp5pre4
    @inbounds @.. broadcast=false k[4]*b40diff4invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                  T::Type{Val{4}})
    @dp5pre4
    @views @.. broadcast=false k[4][idxs]*b40diff4invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs::Nothing,
                                   T::Type{Val{4}})
    @dp5pre4
    @inbounds @.. broadcast=false out=k[4] * b40diff4invdt3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP5ConstantCache, DP5Cache}, idxs,
                                   T::Type{Val{4}})
    @dp5pre4
    @views @.. broadcast=false out=k[4][idxs] * b40diff4invdt3
    out
end

"""
Second order strong stability preserving (SSP) interpolant.

Ketcheson, Lóczi, Jangabylova, Kusmanov: Dense output for SSP RK methods (2017).
"""
@def ssprkpre0 begin
    c00 = @evalpoly(Θ, 1, 0, -1)
    c10 = Θ^2
    b10dt = @evalpoly(Θ, 0, 1, -1) * dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                               SSPRK33ConstantCache, SSPRK33Cache,
                                               SSPRK43ConstantCache, SSPRK43Cache,
                                               SSPRK432ConstantCache, SSPRK432Cache},
                                  idxs::Nothing, T::Type{Val{0}})
    @ssprkpre0
    @inbounds @.. broadcast=false y₀*c00+y₁*c10+k[1]*b10dt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                               SSPRK33ConstantCache, SSPRK33Cache,
                                               SSPRK43ConstantCache, SSPRK43Cache,
                                               SSPRK432ConstantCache, SSPRK432Cache}, idxs,
                                  T::Type{Val{0}})
    @ssprkpre0
    @views @.. broadcast=false y₀[idxs]*c00+y₁[idxs]*c10+k[1][idxs]*b10dt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                                SSPRK33ConstantCache, SSPRK33Cache,
                                                SSPRK43ConstantCache, SSPRK43Cache,
                                                SSPRK432ConstantCache, SSPRK432Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @ssprkpre0
    @inbounds @.. broadcast=false out=y₀ * c00 + y₁ * c10 + k[1] * b10dt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                                SSPRK33ConstantCache, SSPRK33Cache,
                                                SSPRK43ConstantCache, SSPRK43Cache,
                                                SSPRK432ConstantCache, SSPRK432Cache}, idxs,
                                   T::Type{Val{0}})
    @ssprkpre0
    @views @.. broadcast=false out=y₀[idxs] * c00 + y₁[idxs] * c10 + k[1][idxs] * b10dt
    out
end

@def ssprkpre1 begin
    b10diff = @evalpoly(Θ, 1, -2)
    c10diffinvdt = 2Θ * inv(dt) # = -c00diff * inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                               SSPRK33ConstantCache, SSPRK33Cache,
                                               SSPRK43ConstantCache, SSPRK43Cache,
                                               SSPRK432ConstantCache, SSPRK432Cache},
                                  idxs::Nothing, T::Type{Val{1}})
    @ssprkpre1
    @inbounds @.. broadcast=false (y₁ - y₀) * c10diffinvdt+k[1] * b10diff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                               SSPRK33ConstantCache, SSPRK33Cache,
                                               SSPRK43ConstantCache, SSPRK43Cache,
                                               SSPRK432ConstantCache, SSPRK432Cache}, idxs,
                                  T::Type{Val{1}})
    @ssprkpre1
    @views @.. broadcast=false (y₁[idxs] - y₀[idxs]) * c10diffinvdt+k[1][idxs] * b10diff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                                SSPRK33ConstantCache, SSPRK33Cache,
                                                SSPRK43ConstantCache, SSPRK43Cache,
                                                SSPRK432ConstantCache, SSPRK432Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @ssprkpre1
    @inbounds @.. broadcast=false out=(y₁ - y₀) * c10diffinvdt + k[1] * b10diff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                                SSPRK33ConstantCache, SSPRK33Cache,
                                                SSPRK43ConstantCache, SSPRK43Cache,
                                                SSPRK432ConstantCache, SSPRK432Cache}, idxs,
                                   T::Type{Val{1}})
    @ssprkpre1
    @views @.. broadcast=false out=(y₁[idxs] - y₀[idxs]) * c10diffinvdt +
                                   k[1][idxs] * b10diff
    out
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
                                  idxs::Nothing, T::Type{Val{2}})
    @ssprkpre2
    @inbounds @.. broadcast=false (y₁ - y₀) * c10diff2invdt2+k[1] * b10diff2invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                               SSPRK33ConstantCache, SSPRK33Cache,
                                               SSPRK43ConstantCache, SSPRK43Cache,
                                               SSPRK432ConstantCache, SSPRK432Cache}, idxs,
                                  T::Type{Val{2}})
    @ssprkpre2
    @views @.. broadcast=false (y₁[idxs] - y₀[idxs]) *
                               c10diff2invdt2+k[1][idxs] * b10diff2invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                                SSPRK33ConstantCache, SSPRK33Cache,
                                                SSPRK43ConstantCache, SSPRK43Cache,
                                                SSPRK432ConstantCache, SSPRK432Cache},
                                   idxs::Nothing, T::Type{Val{2}})
    @ssprkpre2
    @inbounds @.. broadcast=false out=(y₁ - y₀) * c10diff2invdt2 + k[1] * b10diff2invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{SSPRK22ConstantCache, SSPRK22Cache,
                                                SSPRK33ConstantCache, SSPRK33Cache,
                                                SSPRK43ConstantCache, SSPRK43Cache,
                                                SSPRK432ConstantCache, SSPRK432Cache}, idxs,
                                   T::Type{Val{2}})
    @ssprkpre2
    @views @.. broadcast=false out=(y₁[idxs] - y₀[idxs]) * c10diff2invdt2 +
                                   k[1][idxs] * b10diff2invdt
    out
end

"""
Runge–Kutta pairs of order 5(4) satisfying only the first column
simplifying assumption

Ch. Tsitouras
"""
@def tsit5unpack begin if typeof(cache) <: OrdinaryDiffEqMutableCache
    @unpack r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62, r63, r64, r72, r73, r74 = cache.tab
else
    @unpack r11, r12, r13, r14, r22, r23, r24, r32, r33, r34, r42, r43, r44, r52, r53, r54, r62, r63, r64, r72, r73, r74 = cache
end end

@def tsit5pre0 begin
    @tsit5unpack
    b1Θ = @evalpoly(Θ, 0, r11, r12, r13, r14)
    b2Θ = @evalpoly(Θ, 0, 0, r22, r23, r24)
    b3Θ = @evalpoly(Θ, 0, 0, r32, r33, r34)
    b4Θ = @evalpoly(Θ, 0, 0, r42, r43, r44)
    b5Θ = @evalpoly(Θ, 0, 0, r52, r53, r54)
    b6Θ = @evalpoly(Θ, 0, 0, r62, r63, r64)
    b7Θ = @evalpoly(Θ, 0, 0, r72, r73, r74)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Tsit5ConstantCache,
                                  idxs::Nothing, T::Type{Val{0}})
    @tsit5pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[2]*b2Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ +
                      k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Tsit5Cache, idxs::Nothing,
                                  T::Type{Val{0}})
    @tsit5pre0
    return @inbounds @.. broadcast=false y₀+dt * (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ +
                                             k[4] * b4Θ +
                                             k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                  T::Type{Val{0}})
    @tsit5pre0
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[2][idxs] * b2Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @tsit5pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ +
                                       k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ)
    out
end

@muladd function _ode_interpolant!(out::Array, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @tsit5pre0
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = y₀[i] +
                 dt * (k[1][i] * b1Θ + k[2][i] * b2Θ + k[3][i] * b3Θ + k[4][i] * b4Θ +
                  k[5][i] * b5Θ + k[6][i] * b6Θ + k[7][i] * b7Θ)
    end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                   T::Type{Val{0}})
    @tsit5pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[2][idxs] * b2Θ + k[3][idxs] * b3Θ +
                                    k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ +
                                    k[7][idxs] * b7Θ)
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*(k[1][i]*b1Θ + k[2][i]*b2Θ + k[3][i]*b3Θ + k[4][i]*b4Θ + k[5][i]*b5Θ + k[6][i]*b6Θ + k[7][i]*b7Θ)
    #end
    out
end

@muladd function _ode_interpolant!(out::Array, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                   T::Type{Val{0}})
    @tsit5pre0
    @inbounds for (j, i) in enumerate(idxs)
        out[j] = y₀[i] +
                 dt * (k[1][i] * b1Θ + k[2][i] * b2Θ + k[3][i] * b3Θ + k[4][i] * b4Θ +
                  k[5][i] * b5Θ + k[6][i] * b6Θ + k[7][i] * b7Θ)
    end
    out
end

@def tsit5pre1 begin
    @tsit5unpack
    b1Θdiff = @evalpoly(Θ, r11, 2*r12, 3*r13, 4*r14)
    b2Θdiff = @evalpoly(Θ, 0, 2*r22, 3*r23, 4*r24)
    b3Θdiff = @evalpoly(Θ, 0, 2*r32, 3*r33, 4*r34)
    b4Θdiff = @evalpoly(Θ, 0, 2*r42, 3*r43, 4*r44)
    b5Θdiff = @evalpoly(Θ, 0, 2*r52, 3*r53, 4*r54)
    b6Θdiff = @evalpoly(Θ, 0, 2*r62, 3*r63, 4*r64)
    b7Θdiff = @evalpoly(Θ, 0, 2*r72, 3*r73, 4*r74)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Tsit5ConstantCache,
                                  idxs::Nothing, T::Type{Val{1}})
    @tsit5pre1
    # return @.. broadcast=false k[1]*b1Θdiff + k[2]*b2Θdiff + k[3]*b3Θdiff + k[4]*b4Θdiff + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff
    return @inbounds k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
                     k[5] * b5Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Tsit5Cache, idxs::Nothing,
                                  T::Type{Val{1}})
    @tsit5pre1
    return @inbounds @.. broadcast=false k[1]*b1Θdiff+k[2]*b2Θdiff+k[3]*b3Θdiff+
                                         k[4]*b4Θdiff+k[5]*b5Θdiff+k[6]*b6Θdiff+k[7]*b7Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                  T::Type{Val{1}})
    @tsit5pre1
    # return @.. broadcast=false k[1][idxs]*b1Θdiff + k[2][idxs]*b2Θdiff + k[3][idxs]*b3Θdiff + k[4][idxs]*b4Θdiff + k[5][idxs]*b5Θdiff + k[6][idxs]*b6Θdiff + k[7][idxs]*b7Θdiff
    return k[1][idxs] * b1Θdiff + k[2][idxs] * b2Θdiff + k[3][idxs] * b3Θdiff +
           k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
           k[7][idxs] * b7Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @tsit5pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff +
                                      k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
                                      k[7] * b7Θdiff
    #@inbounds for i in eachindex(out)
    #  out[i] = k[1][i]*b1Θdiff + k[2][i]*b2Θdiff + k[3][i]*b3Θdiff + k[4][i]*b4Θdiff + k[5][i]*b5Θdiff + k[6][i]*b6Θdiff + k[7][i]*b7Θdiff
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                   T::Type{Val{1}})
    @tsit5pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[2][idxs] * b2Θdiff +
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
    b1Θdiff2 = @evalpoly(Θ, 2*r12, 6*r13, 12*r14)
    b2Θdiff2 = @evalpoly(Θ, 2*r22, 6*r23, 12*r24)
    b3Θdiff2 = @evalpoly(Θ, 2*r32, 6*r33, 12*r34)
    b4Θdiff2 = @evalpoly(Θ, 2*r42, 6*r43, 12*r44)
    b5Θdiff2 = @evalpoly(Θ, 2*r52, 6*r53, 12*r54)
    b6Θdiff2 = @evalpoly(Θ, 2*r62, 6*r63, 12*r64)
    b7Θdiff2 = @evalpoly(Θ, 2*r72, 6*r73, 12*r74)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                  idxs::Nothing, T::Type{Val{2}})
    @tsit5pre2
    # return @.. broadcast=false k[1]*b1Θdiff2 + k[2]*b2Θdiff2 + k[3]*b3Θdiff2 + k[4]*b4Θdiff2 + k[5]*b5Θdiff2 + k[6]*b6Θdiff2 + k[7]*b7Θdiff2
    return @inbounds (k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
                      k[4] * b4Θdiff2 +
                      k[5] * b5Θdiff2 + k[6] * b6Θdiff2 + k[7] * b7Θdiff2) * invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                  T::Type{Val{2}})
    @tsit5pre2
    # return @.. broadcast=false k[1][idxs]*b1Θdiff2 + k[2][idxs]*b2Θdiff2 + k[3][idxs]*b3Θdiff2 + k[4][idxs]*b4Θdiff2 + k[5][idxs]*b5Θdiff2 + k[6][idxs]*b6Θdiff2 + k[7][idxs]*b7Θdiff2
    return (k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 + k[3][idxs] * b3Θdiff2 +
            k[4][idxs] * b4Θdiff2 + k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2 +
            k[7][idxs] * b7Θdiff2) * invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                   idxs::Nothing, T::Type{Val{2}})
    @tsit5pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
                                       k[4] * b4Θdiff2 + k[5] * b5Θdiff2 + k[6] * b6Θdiff2 +
                                       k[7] * b7Θdiff2) * invdt
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff2 + k[2][i]*b2Θdiff2 + k[3][i]*b3Θdiff2 + k[4][i]*b4Θdiff2 + k[5][i]*b5Θdiff2 + k[6][i]*b6Θdiff2 + k[7][i]*b7Θdiff2)*invdt
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                   T::Type{Val{2}})
    @tsit5pre2
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 +
                                    k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2 +
                                    k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2 +
                                    k[7][idxs] * b7Θdiff2) * invdt
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff2 + k[2][i]*b2Θdiff2 + k[3][i]*b3Θdiff2 + k[4][i]*b4Θdiff2 + k[5][i]*b5Θdiff2 + k[6][i]*b6Θdiff2 + k[7][i]*b7Θdiff2)*invdt
    #end
    out
end

@def tsit5pre3 begin
    @tsit5unpack
    b1Θdiff3 = @evalpoly(Θ, 6*r13, 24*r14)
    b2Θdiff3 = @evalpoly(Θ, 6*r23, 24*r24)
    b3Θdiff3 = @evalpoly(Θ, 6*r33, 24*r34)
    b4Θdiff3 = @evalpoly(Θ, 6*r43, 24*r44)
    b5Θdiff3 = @evalpoly(Θ, 6*r53, 24*r54)
    b6Θdiff3 = @evalpoly(Θ, 6*r63, 24*r64)
    b7Θdiff3 = @evalpoly(Θ, 6*r73, 24*r74)
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                  idxs::Nothing, T::Type{Val{3}})
    @tsit5pre3
    # return @.. broadcast=false k[1]*b1Θdiff3 + k[2]*b2Θdiff3 + k[3]*b3Θdiff3 + k[4]*b4Θdiff3 + k[5]*b5Θdiff3 + k[6]*b6Θdiff3 + k[7]*b7Θdiff3
    return @inbounds (k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
                      k[4] * b4Θdiff3 +
                      k[5] * b5Θdiff3 + k[6] * b6Θdiff3 + k[7] * b7Θdiff3) * invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                  T::Type{Val{3}})
    @tsit5pre3
    # return @.. broadcast=false k[1][idxs]*b1Θdiff3 + k[2][idxs]*b2Θdiff3 + k[3][idxs]*b3Θdiff3 + k[4][idxs]*b4Θdiff3 + k[5][idxs]*b5Θdiff3 + k[6][idxs]*b6Θdiff3 + k[7][idxs]*b7Θdiff3
    return (k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 + k[3][idxs] * b3Θdiff3 +
            k[4][idxs] * b4Θdiff3 + k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3 +
            k[7][idxs] * b7Θdiff3) * invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                   idxs::Nothing, T::Type{Val{3}})
    @tsit5pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
                                       k[4] * b4Θdiff3 + k[5] * b5Θdiff3 + k[6] * b6Θdiff3 +
                                       k[7] * b7Θdiff3) * invdt2
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff3 + k[2][i]*b2Θdiff3 + k[3][i]*b3Θdiff3 + k[4][i]*b4Θdiff3 + k[5][i]*b5Θdiff3 + k[6][i]*b6Θdiff3 + k[7][i]*b7Θdiff3)*invdt2
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                   T::Type{Val{3}})
    @tsit5pre3
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 +
                                    k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3 +
                                    k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3 +
                                    k[7][idxs] * b7Θdiff3) * invdt2
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

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                  idxs::Nothing, T::Type{Val{4}})
    @tsit5pre4
    # return @.. broadcast=false k[1]*b1Θdiff4 + k[2]*b2Θdiff4 + k[3]*b3Θdiff4 + k[4]*b4Θdiff4 + k[5]*b5Θdiff4 + k[6]*b6Θdiff4 + k[7]*b7Θdiff4
    return @inbounds (k[1] * b1Θdiff4 + k[2] * b2Θdiff4 + k[3] * b3Θdiff4 +
                      k[4] * b4Θdiff4 +
                      k[5] * b5Θdiff4 + k[6] * b6Θdiff4 + k[7] * b7Θdiff4) * invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                  T::Type{Val{4}})
    @tsit5pre4
    # return @.. broadcast=false k[1][idxs]*b1Θdiff4 + k[2][idxs]*b2Θdiff4 + k[3][idxs]*b3Θdiff4 + k[4][idxs]*b4Θdiff4 + k[5][idxs]*b5Θdiff4 + k[6][idxs]*b6Θdiff4 + k[7][idxs]*b7Θdiff4
    return (k[1][idxs] * b1Θdiff4 + k[2][idxs] * b2Θdiff4 + k[3][idxs] * b3Θdiff4 +
            k[4][idxs] * b4Θdiff4 + k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4 +
            k[7][idxs] * b7Θdiff4) * invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache},
                                   idxs::Nothing, T::Type{Val{4}})
    @tsit5pre4
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff4 + k[2] * b2Θdiff4 + k[3] * b3Θdiff4 +
                                       k[4] * b4Θdiff4 + k[5] * b5Θdiff4 + k[6] * b6Θdiff4 +
                                       k[7] * b7Θdiff4) * invdt3
    #@inbounds for i in eachindex(out)
    #  out[i] = (k[1][i]*b1Θdiff4 + k[2][i]*b2Θdiff4 + k[3][i]*b3Θdiff4 + k[4][i]*b4Θdiff4 + k[5][i]*b5Θdiff4 + k[6][i]*b6Θdiff4 + k[7][i]*b7Θdiff4)*invdt3
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Tsit5ConstantCache, Tsit5Cache}, idxs,
                                   T::Type{Val{4}})
    @tsit5pre4
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff4 + k[2][idxs] * b2Θdiff4 +
                                    k[3][idxs] * b3Θdiff4 + k[4][idxs] * b4Θdiff4 +
                                    k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4 +
                                    k[7][idxs] * b7Θdiff4) * invdt3
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = (k[1][i]*b1Θdiff4 + k[2][i]*b2Θdiff4 + k[3][i]*b3Θdiff4 + k[4][i]*b4Θdiff4 + k[5][i]*b5Θdiff4 + k[6][i]*b6Θdiff4 + k[7][i]*b7Θdiff4)*invdt3
    #end
    out
end

"""

"""
@def owrenzen3unpack begin if typeof(cache) <: OrdinaryDiffEqMutableCache
    @unpack r13, r12, r23, r22, r33, r32 = cache.tab
else
    @unpack r13, r12, r23, r22, r33, r32 = cache
end end

@def owrenzen3pre0 begin
    @owrenzen3unpack
    b1Θ = @evalpoly(Θ, 0, 1, r12, r13)
    b2Θ = @evalpoly(Θ, 0, 0, r22, r23)
    b3Θ = @evalpoly(Θ, 0, 0, r32, r33)
    b4Θ = @evalpoly(Θ, 0, 0, -1, 1)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                  idxs::Nothing, T::Type{Val{0}})
    @owrenzen3pre0
    @inbounds @.. broadcast=false y₀+dt *
                                     (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                  idxs, T::Type{Val{0}})
    @owrenzen3pre0
    @views @.. broadcast=false y₀[idxs]+dt * (k[1][idxs] * b1Θ + k[2][idxs] * b2Θ +
                                         k[3][idxs] * b3Θ +
                                         k[4][idxs] * b4Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @owrenzen3pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[2] * b2Θ + k[3] * b3Θ + k[4] * b4Θ)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs, T::Type{Val{0}})
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
    b2Θdiff = @evalpoly(Θ, 0, 2*r22, 3*r23)
    b3Θdiff = @evalpoly(Θ, 0, 2*r32, 3*r33)
    b4Θdiff = @evalpoly(Θ, 0, -2, 3)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                  idxs::Nothing, T::Type{Val{1}})
    @owrenzen3pre1
    @inbounds @.. broadcast=false k[1]*b1Θdiff+k[2]*b2Θdiff+k[3]*b3Θdiff+k[4]*b4Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                  idxs, T::Type{Val{1}})
    @owrenzen3pre1
    @views @.. broadcast=false k[1][idxs]*b1Θdiff+k[2][idxs]*b2Θdiff+k[3][idxs]*b3Θdiff+
                               k[4][idxs]*b4Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @owrenzen3pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[2] * b2Θdiff + k[3] * b3Θdiff +
                                      k[4] * b4Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs, T::Type{Val{1}})
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
                                  idxs::Nothing, T::Type{Val{2}})
    @owrenzen3pre2
    @inbounds @.. broadcast=false (k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
                                   k[4] * b4Θdiff2)*invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                  idxs, T::Type{Val{2}})
    @owrenzen3pre2
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff2 + k[2][idxs] * b2Θdiff2 +
                                k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2)*invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs::Nothing, T::Type{Val{2}})
    @owrenzen3pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff2 + k[2] * b2Θdiff2 + k[3] * b3Θdiff2 +
                                       k[4] * b4Θdiff2) * invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs, T::Type{Val{2}})
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
                                  idxs::Nothing, T::Type{Val{3}})
    @owrenzen3pre3
    @inbounds @.. broadcast=false (k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
                                   k[4] * b4Θdiff3)*invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                  idxs, T::Type{Val{3}})
    @owrenzen3pre3
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 +
                                k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3)*invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs::Nothing, T::Type{Val{3}})
    @owrenzen3pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff3 + k[2] * b2Θdiff3 + k[3] * b3Θdiff3 +
                                       k[4] * b4Θdiff3) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen3ConstantCache, OwrenZen3Cache},
                                   idxs, T::Type{Val{3}})
    @owrenzen3pre3
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff3 + k[2][idxs] * b2Θdiff3 +
                                    k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3) * invdt2
    out
end

"""

"""
@def owrenzen4unpack begin if typeof(cache) <: OrdinaryDiffEqMutableCache
    @unpack r14, r13, r12, r34, r33, r32, r44, r43, r42, r54, r53, r52, r64, r63, r62 = cache.tab
else
    @unpack r14, r13, r12, r34, r33, r32, r44, r43, r42, r54, r53, r52, r64, r63, r62 = cache
end end

@def owrenzen4pre0 begin
    @owrenzen4unpack
    b1Θ = @evalpoly(Θ, 0, 1, r12, r13, r14)
    b3Θ = @evalpoly(Θ, 0, 0, r32, r33, r34)
    b4Θ = @evalpoly(Θ, 0, 0, r42, r43, r44)
    b5Θ = @evalpoly(Θ, 0, 0, r52, r53, r54)
    b6Θ = @evalpoly(Θ, 0, 0, r62, r63, r64)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs::Nothing, T::Type{Val{0}})
    @owrenzen4pre0
    # return @.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs, T::Type{Val{0}})
    @owrenzen4pre0
    # return @.. broadcast=false y₀[idxs] + dt*(k[1][idxs]*b1Θ + k[3][idxs]*b3Θ +
    #                          k[4][idxs]*b4Θ + k[5][idxs]*b5Θ + k[6][idxs]*b6Θ)
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[3][idxs] * b3Θ +
            k[4][idxs] * b4Θ + k[5][idxs] * b5Θ + k[6][idxs] * b6Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs::Nothing, T::Type{Val{0}})
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
                                   idxs, T::Type{Val{0}})
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
    b3Θdiff = @evalpoly(Θ, 0, 2*r32, 3*r33, 4*r34)
    b4Θdiff = @evalpoly(Θ, 0, 2*r42, 3*r43, 4*r44)
    b5Θdiff = @evalpoly(Θ, 0, 2*r52, 3*r53, 4*r54)
    b6Θdiff = @evalpoly(Θ, 0, 2*r62, 3*r63, 4*r64)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs::Nothing, T::Type{Val{1}})
    @owrenzen4pre1
    @inbounds @.. broadcast=false k[1]*b1Θdiff+k[3]*b3Θdiff+k[4]*b4Θdiff+k[5]*b5Θdiff+
                                  k[6]*b6Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs, T::Type{Val{1}})
    @owrenzen4pre1
    @views @.. broadcast=false k[1][idxs]*b1Θdiff+k[3][idxs]*b3Θdiff+k[4][idxs]*b4Θdiff+
                               k[5][idxs]*b5Θdiff+k[6][idxs]*b6Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @owrenzen4pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
                                      k[5] * b5Θdiff + k[6] * b6Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs, T::Type{Val{1}})
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
                                  idxs::Nothing, T::Type{Val{2}})
    @owrenzen4pre2
    @.. broadcast=false (k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                         k[5] * b5Θdiff2 + k[6] * b6Θdiff2)*invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs, T::Type{Val{2}})
    @owrenzen4pre2
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff2 + k[3][idxs] * b3Θdiff2 +
                                k[4][idxs] * b4Θdiff2 +
                                k[5][idxs] * b5Θdiff2 + k[6][idxs] * b6Θdiff2)*invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs::Nothing, T::Type{Val{2}})
    @owrenzen4pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                                       k[5] * b5Θdiff2 + k[6] * b6Θdiff2) * invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs, T::Type{Val{2}})
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
                                  idxs::Nothing, T::Type{Val{3}})
    @owrenzen4pre3
    @inbounds @.. broadcast=false (k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                                   k[5] * b5Θdiff3 + k[6] * b6Θdiff3)*invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs, T::Type{Val{3}})
    @owrenzen4pre3
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff3 + k[3][idxs] * b3Θdiff3 +
                                k[4][idxs] * b4Θdiff3 +
                                k[5][idxs] * b5Θdiff3 + k[6][idxs] * b6Θdiff3)*invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs::Nothing, T::Type{Val{3}})
    @owrenzen4pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                                       k[5] * b5Θdiff3 + k[6] * b6Θdiff3) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs, T::Type{Val{3}})
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
                                  idxs::Nothing, T::Type{Val{4}})
    @owrenzen4pre4
    @.. broadcast=false (k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                         k[5] * b5Θdiff4 + k[6] * b6Θdiff4)*invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                  idxs, T::Type{Val{4}})
    @owrenzen4pre4
    @views @.. broadcast=false (k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 +
                                k[4][idxs] * b4Θdiff4 +
                                k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4)*invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs::Nothing, T::Type{Val{4}})
    @owrenzen4pre4
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                                       k[5] * b5Θdiff4 + k[6] * b6Θdiff4) * invdt3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen4ConstantCache, OwrenZen4Cache},
                                   idxs, T::Type{Val{4}})
    @owrenzen4pre4
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 +
                                    k[4][idxs] * b4Θdiff4 +
                                    k[5][idxs] * b5Θdiff4 + k[6][idxs] * b6Θdiff4) * invdt3
    out
end

"""

"""
@def owrenzen5unpack begin if typeof(cache) <: OrdinaryDiffEqMutableCache
    @unpack r15, r14, r13, r12, r35, r34, r33, r32, r45, r44, r43, r42, r55, r54, r53, r52, r65, r64, r63, r62, r75, r74, r73, r72, r85, r84, r83, r82 = cache.tab
else
    @unpack r15, r14, r13, r12, r35, r34, r33, r32, r45, r44, r43, r42, r55, r54, r53, r52, r65, r64, r63, r62, r75, r74, r73, r72, r85, r84, r83, r82 = cache
end end

@def owrenzen5pre0 begin
    @owrenzen5unpack
    b1Θ = @evalpoly(Θ, 0, 1, r12, r13, r14, r15)
    b3Θ = @evalpoly(Θ, 0, 0, r32, r33, r34, r35)
    b4Θ = @evalpoly(Θ, 0, 0, r42, r43, r44, r45)
    b5Θ = @evalpoly(Θ, 0, 0, r52, r53, r54, r55)
    b6Θ = @evalpoly(Θ, 0, 0, r62, r63, r64, r65)
    b7Θ = @evalpoly(Θ, 0, 0, r72, r73, r74, r75)
    b8Θ = @evalpoly(Θ, 0, 0, r82, r83, r84, r85)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs::Nothing, T::Type{Val{0}})
    @owrenzen5pre0
    # return @.. broadcast=false y₀ + dt*(k[1]*b1Θ  + k[3]*b3Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ +
    #                    k[7]*b7Θ + k[8]*b8Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ +
                      k[7] * b7Θ + k[8] * b8Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs, T::Type{Val{0}})
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
                                   idxs::Nothing, T::Type{Val{0}})
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
                                   idxs, T::Type{Val{0}})
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
    b3Θdiff = @evalpoly(Θ, 0, 2*r32, 3*r33, 4*r34, 5*r35)
    b4Θdiff = @evalpoly(Θ, 0, 2*r42, 3*r43, 4*r44, 5*r45)
    b5Θdiff = @evalpoly(Θ, 0, 2*r52, 3*r53, 4*r54, 5*r55)
    b6Θdiff = @evalpoly(Θ, 0, 2*r62, 3*r63, 4*r64, 5*r65)
    b7Θdiff = @evalpoly(Θ, 0, 2*r72, 3*r73, 4*r74, 5*r75)
    b8Θdiff = @evalpoly(Θ, 0, 2*r82, 3*r83, 4*r84, 5*r85)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs::Nothing, T::Type{Val{1}})
    @owrenzen5pre1
    return @inbounds k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff +
                     k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs, T::Type{Val{1}})
    @owrenzen5pre1
    k[1][idxs] * b1Θdiff + k[3][idxs] * b3Θdiff + k[4][idxs] * b4Θdiff +
    k[5][idxs] * b5Θdiff +
    k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                   idxs::Nothing, T::Type{Val{1}})
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
                                   idxs, T::Type{Val{1}})
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
                                  idxs::Nothing, T::Type{Val{2}})
    @owrenzen5pre2
    return @inbounds (k[1] * b1Θdiff2 + k[3] * b3Θdiff2 + k[4] * b4Θdiff2 +
                      k[5] * b5Θdiff2 +
                      k[6] * b6Θdiff2 + k[7] * b7Θdiff2 + k[8] * b8Θdiff2) * invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs, T::Type{Val{2}})
    @owrenzen5pre2
    (k[1][idxs] * b1Θdiff2 + k[3][idxs] * b3Θdiff2 + k[4][idxs] * b4Θdiff2 +
     k[5][idxs] * b5Θdiff2 +
     k[6][idxs] * b6Θdiff2 + k[7][idxs] * b7Θdiff2 + k[8][idxs] * b8Θdiff2) * invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                   idxs::Nothing, T::Type{Val{2}})
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
                                   idxs, T::Type{Val{2}})
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
                                  idxs::Nothing, T::Type{Val{3}})
    @owrenzen5pre3
    return @inbounds (k[1] * b1Θdiff3 + k[3] * b3Θdiff3 + k[4] * b4Θdiff3 +
                      k[5] * b5Θdiff3 +
                      k[6] * b6Θdiff3 + k[7] * b7Θdiff3 + k[8] * b8Θdiff3) * invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs, T::Type{Val{3}})
    @owrenzen5pre3
    (k[1][idxs] * b1Θdiff3 + k[3][idxs] * b3Θdiff3 + k[4][idxs] * b4Θdiff3 +
     k[5][idxs] * b5Θdiff3 +
     k[6][idxs] * b6Θdiff3 + k[7][idxs] * b7Θdiff3 + k[8][idxs] * b8Θdiff3) * invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                   idxs::Nothing, T::Type{Val{3}})
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
                                   idxs, T::Type{Val{3}})
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
                                  idxs::Nothing, T::Type{Val{4}})
    @owrenzen5pre4
    return @inbounds (k[1] * b1Θdiff4 + k[3] * b3Θdiff4 + k[4] * b4Θdiff4 +
                      k[5] * b5Θdiff4 +
                      k[6] * b6Θdiff4 + k[7] * b7Θdiff4 + k[8] * b8Θdiff4) * invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs, T::Type{Val{4}})
    @owrenzen5pre4
    (k[1][idxs] * b1Θdiff4 + k[3][idxs] * b3Θdiff4 + k[4][idxs] * b4Θdiff4 +
     k[5][idxs] * b5Θdiff4 +
     k[6][idxs] * b6Θdiff4 + k[7][idxs] * b7Θdiff4 + k[8][idxs] * b8Θdiff4) * invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                   idxs::Nothing, T::Type{Val{4}})
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
                                   idxs, T::Type{Val{4}})
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
                                  idxs::Nothing, T::Type{Val{5}})
    @owrenzen5pre5
    return @inbounds (k[1] * b1Θdiff5 + k[3] * b3Θdiff5 + k[4] * b4Θdiff5 +
                      k[5] * b5Θdiff5 +
                      k[6] * b6Θdiff5 + k[7] * b7Θdiff5 + k[8] * b8Θdiff5) * invdt4
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                  idxs, T::Type{Val{5}})
    @owrenzen5pre5
    (k[1][idxs] * b1Θdiff5 + k[3][idxs] * b3Θdiff5 + k[4][idxs] * b4Θdiff5 +
     k[5][idxs] * b5Θdiff5 +
     k[6][idxs] * b6Θdiff5 + k[7][idxs] * b7Θdiff5 + k[8][idxs] * b8Θdiff5) * invdt4
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{OwrenZen5ConstantCache, OwrenZen5Cache},
                                   idxs::Nothing, T::Type{Val{5}})
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
                                   idxs, T::Type{Val{5}})
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
@def bs5unpack begin if typeof(cache) <: OrdinaryDiffEqMutableCache
    @unpack r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115, r114, r113, r112 = cache.tab
else
    @unpack r016, r015, r014, r013, r012, r036, r035, r034, r033, r032, r046, r045, r044, r043, r042, r056, r055, r054, r053, r052, r066, r065, r064, r063, r062, r076, r075, r074, r073, r072, r086, r085, r084, r083, r082, r096, r095, r094, r093, r106, r105, r104, r103, r102, r116, r115, r114, r113, r112 = cache
end end

@def bs5pre0 begin
    @bs5unpack
    b1Θ = @evalpoly(Θ, 0, 0, r012, r013, r014, r015, r016)
    b3Θ = @evalpoly(Θ, 0, 0, r032, r033, r034, r035, r036)
    b4Θ = @evalpoly(Θ, 0, 0, r042, r043, r044, r045, r046)
    b5Θ = @evalpoly(Θ, 0, 0, r052, r053, r054, r055, r056)
    b6Θ = @evalpoly(Θ, 0, 0, r062, r063, r064, r065, r066)
    b7Θ = @evalpoly(Θ, 0, 0, r072, r073, r074, r075, r076)
    b8Θ = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086)
    b9Θ = @evalpoly(Θ, 0, 0, 0, r093, r094, r095, r096)
    b10Θ = @evalpoly(Θ, 0, 0, r102, r103, r104, r105, r106)
    b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::BS5ConstantCache, idxs::Nothing,
                                  T::Type{Val{0}})
    @bs5pre0
    # return @.. broadcast=false y₀ + dt*Θ*k[1] + dt*(k[1]*b1Θ  + k[3]*b3Θ + k[4]*b4Θ  + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ)
    return @inbounds y₀ + dt * Θ * k[1] +
                     dt * (k[1] * b1Θ + k[3] * b3Θ + k[4] * b4Θ + k[5] * b5Θ +
                      k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ +
                      k[11] * b11Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::BS5Cache, idxs::Nothing,
                                  T::Type{Val{0}})
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
                                  T::Type{Val{0}})
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
                                   T::Type{Val{0}})
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
                                   T::Type{Val{0}})
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
    b1Θdiff = @evalpoly(Θ, 0, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016)
    b3Θdiff = @evalpoly(Θ, 0, 2*r032, 3*r033, 4*r034, 5*r035, 6*r036)
    b4Θdiff = @evalpoly(Θ, 0, 2*r042, 3*r043, 4*r044, 5*r045, 6*r046)
    b5Θdiff = @evalpoly(Θ, 0, 2*r052, 3*r053, 4*r054, 5*r055, 6*r056)
    b6Θdiff = @evalpoly(Θ, 0, 2*r062, 3*r063, 4*r064, 5*r065, 6*r066)
    b7Θdiff = @evalpoly(Θ, 0, 2*r072, 3*r073, 4*r074, 5*r075, 6*r076)
    b8Θdiff = @evalpoly(Θ, 0, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086)
    b9Θdiff = @evalpoly(Θ, 0, 0, 3*r093, 4*r094, 5*r095, 6*r096)
    b10Θdiff = @evalpoly(Θ, 0, 2*r102, 3*r103, 4*r104, 5*r105, 6*r106)
    b11Θdiff = @evalpoly(Θ, 0, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{BS5ConstantCache, BS5Cache}, idxs::Nothing,
                                  T::Type{Val{1}})
    @bs5pre1
    # return @.. broadcast=false k[1] + k[1]*b1Θdiff  + k[3]*b3Θdiff + k[4]*b4Θdiff  + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[10]*b10Θdiff + k[11]*b11Θdiff
    return @inbounds k[1] + k[1] * b1Θdiff + k[3] * b3Θdiff + k[4] * b4Θdiff +
                     k[5] * b5Θdiff +
                     k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff +
                     k[10] * b10Θdiff + k[11] * b11Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{BS5ConstantCache, BS5Cache}, idxs,
                                  T::Type{Val{1}})
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
                                   T::Type{Val{1}})
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
                                   T::Type{Val{1}})
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

## Vern6
@def vern6unpack begin
    @unpack r011, r012, r013, r014, r015, r016, r042, r043, r044, r045, r046, r052, r053, r054, r055, r056, r062, r063, r064, r065, r066, r072, r073, r074, r075, r076, r082, r083, r084, r085, r086, r092, r093, r094, r095, r096, r102, r103, r104, r105, r106, r112, r113, r114, r115, r116, r122, r123, r124, r125, r126 = cache.tab.interp
end

@def vern6pre0 begin
    @vern6unpack
    b1Θ = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016)
    b4Θ = @evalpoly(Θ, 0, 0, r042, r043, r044, r045, r046)
    b5Θ = @evalpoly(Θ, 0, 0, r052, r053, r054, r055, r056)
    b6Θ = @evalpoly(Θ, 0, 0, r062, r063, r064, r065, r066)
    b7Θ = @evalpoly(Θ, 0, 0, r072, r073, r074, r075, r076)
    b8Θ = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086)
    b9Θ = @evalpoly(Θ, 0, 0, r092, r093, r094, r095, r096)
    b10Θ = @evalpoly(Θ, 0, 0, r102, r103, r104, r105, r106)
    b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116)
    b12Θ = @evalpoly(Θ, 0, 0, r122, r123, r124, r125, r126)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern6ConstantCache,
                                  idxs::Nothing, T::Type{Val{0}})
    @vern6pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ +
                      k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ + k[11] * b11Θ + k[12] * b12Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern6Cache, idxs::Nothing,
                                  T::Type{Val{0}})
    @vern6pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ)
    return @inbounds @.. broadcast=false y₀+dt * (k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ +
                                             k[6] * b6Θ + k[7] * b7Θ +
                                             k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ +
                                             k[11] * b11Θ + k[12] * b12Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
                                  T::Type{Val{0}})
    @vern6pre0
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
            k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
            k[10][idxs] * b10Θ + k[11][idxs] * b11Θ + k[12][idxs] * b12Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern6ConstantCache, Vern6Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @vern6pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ +
                                       k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ + k[10] * b10Θ +
                                       k[11] * b11Θ + k[12] * b12Θ)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
                                   T::Type{Val{0}})
    @vern6pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
                                    k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ +
                                    k[9][idxs] * b9Θ + k[10][idxs] * b10Θ +
                                    k[11][idxs] * b11Θ + k[12][idxs] * b12Θ)
    out
end

@def vern6pre1 begin
    @vern6unpack
    b1Θdiff = @evalpoly(Θ, r011, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016)
    b4Θdiff = @evalpoly(Θ, 0, 2*r042, 3*r043, 4*r044, 5*r045, 6*r046)
    b5Θdiff = @evalpoly(Θ, 0, 2*r052, 3*r053, 4*r054, 5*r055, 6*r056)
    b6Θdiff = @evalpoly(Θ, 0, 2*r062, 3*r063, 4*r064, 5*r065, 6*r066)
    b7Θdiff = @evalpoly(Θ, 0, 2*r072, 3*r073, 4*r074, 5*r075, 6*r076)
    b8Θdiff = @evalpoly(Θ, 0, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086)
    b9Θdiff = @evalpoly(Θ, 0, 2*r092, 3*r093, 4*r094, 5*r095, 6*r096)
    b10Θdiff = @evalpoly(Θ, 0, 2*r102, 3*r103, 4*r104, 5*r105, 6*r106)
    b11Θdiff = @evalpoly(Θ, 0, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116)
    b12Θdiff = @evalpoly(Θ, 0, 2*r122, 3*r123, 4*r124, 5*r125, 6*r126)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern6ConstantCache, Vern6Cache},
                                  idxs::Nothing, T::Type{Val{1}})
    @vern6pre1
    #@.. broadcast=false k[1]*b1Θdiff + k[4]*b4Θdiff + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[10]*b10Θdiff + k[11]*b11Θdiff + k[12]*b12Θdiff
    return @inbounds k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
                     k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff + k[10] * b10Θdiff +
                     k[11] * b11Θdiff + k[12] * b12Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
                                  T::Type{Val{1}})
    @vern6pre1
    return k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff +
           k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
           k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff + k[11][idxs] * b11Θdiff +
           k[12][idxs] * b12Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern6ConstantCache, Vern6Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @vern6pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff +
                                      k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff +
                                      k[9] * b9Θdiff + k[10] * b10Θdiff + k[11] * b11Θdiff +
                                      k[12] * b12Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern6ConstantCache, Vern6Cache}, idxs,
                                   T::Type{Val{1}})
    @vern6pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff +
                                   k[5][idxs] * b5Θdiff + k[6][idxs] * b6Θdiff +
                                   k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
                                   k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff +
                                   k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff
    out
end

## Vern7
@def vern7unpack begin
    T = Float64
    r011 = convert(T, 1)
    r012 = convert(T, -8.413387198332767)
    r013 = convert(T, 33.675508884490895)
    r014 = convert(T, -70.80159089484886)
    r015 = convert(T, 80.64695108301298)
    r016 = convert(T, -47.19413969837522)
    r017 = convert(T, 11.133813442539243)
    r042 = convert(T, 8.754921980674396)
    r043 = convert(T, -88.4596828699771)
    r044 = convert(T, 346.9017638429916)
    r045 = convert(T, -629.2580030059837)
    r046 = convert(T, 529.6773755604193)
    r047 = convert(T, -167.35886986514018)
    r052 = convert(T, 8.913387586637922)
    r053 = convert(T, -90.06081846893218)
    r054 = convert(T, 353.1807459217058)
    r055 = convert(T, -640.6476819744374)
    r056 = convert(T, 539.2646279047156)
    r057 = convert(T, -170.38809442991547)
    r062 = convert(T, 5.1733120298478)
    r063 = convert(T, -52.271115900055385)
    r064 = convert(T, 204.9853867374073)
    r065 = convert(T, -371.8306118563603)
    r066 = convert(T, 312.9880934374529)
    r067 = convert(T, -98.89290352172495)
    r072 = convert(T, 16.79537744079696)
    r073 = convert(T, -169.70040000059728)
    r074 = convert(T, 665.4937727009246)
    r075 = convert(T, -1207.1638892336007)
    r076 = convert(T, 1016.1291515818546)
    r077 = convert(T, -321.06001557237494)
    r082 = convert(T, -10.005997536098665)
    r083 = convert(T, 101.1005433052275)
    r084 = convert(T, -396.47391512378437)
    r085 = convert(T, 719.1787707014183)
    r086 = convert(T, -605.3681033918824)
    r087 = convert(T, 191.27439892797935)
    r092 = convert(T, 2.764708833638599)
    r093 = convert(T, -27.934602637390462)
    r094 = convert(T, 109.54779186137893)
    r095 = convert(T, -198.7128113064482)
    r096 = convert(T, 167.26633571640318)
    r097 = convert(T, -52.85010499525706)
    r112 = convert(T, -2.1696320280163506)
    r113 = convert(T, 22.016696037569876)
    r114 = convert(T, -86.90152427798948)
    r115 = convert(T, 159.22388973861476)
    r116 = convert(T, -135.9618306534588)
    r117 = convert(T, 43.792401183280006)
    r122 = convert(T, -4.890070188793804)
    r123 = convert(T, 22.75407737425176)
    r124 = convert(T, -30.78034218537731)
    r125 = convert(T, -2.797194317207249)
    r126 = convert(T, 31.369456637508403)
    r127 = convert(T, -15.655927320381801)
    r132 = convert(T, 10.862170929551967)
    r133 = convert(T, -50.542971417827104)
    r134 = convert(T, 68.37148040407511)
    r135 = convert(T, 6.213326521632409)
    r136 = convert(T, -69.68006323194157)
    r137 = convert(T, 34.776056794509195)
    r142 = convert(T, -11.37286691922923)
    r143 = convert(T, 130.79058078246717)
    r144 = convert(T, -488.65113677785604)
    r145 = convert(T, 832.2148793276441)
    r146 = convert(T, -664.7743368554426)
    r147 = convert(T, 201.79288044241662)
    r152 = convert(T, -5.919778732715007)
    r153 = convert(T, 63.27679965889219)
    r154 = convert(T, -265.432682088738)
    r155 = convert(T, 520.1009254140611)
    r156 = convert(T, -467.412109533902)
    r157 = convert(T, 155.3868452824017)
    r162 = convert(T, -10.492146197961823)
    r163 = convert(T, 105.35538525188011)
    r164 = convert(T, -409.43975011988937)
    r165 = convert(T, 732.831448907654)
    r166 = convert(T, -606.3044574733512)
    r167 = convert(T, 188.0495196316683)
end

@def vern7pre0 begin
    @vern7unpack
    b1Θ = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017)
    b4Θ = @evalpoly(Θ, 0, 0, r042, r043, r044, r045, r046, r047)
    b5Θ = @evalpoly(Θ, 0, 0, r052, r053, r054, r055, r056, r057)
    b6Θ = @evalpoly(Θ, 0, 0, r062, r063, r064, r065, r066, r067)
    b7Θ = @evalpoly(Θ, 0, 0, r072, r073, r074, r075, r076, r077)
    b8Θ = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086, r087)
    b9Θ = @evalpoly(Θ, 0, 0, r092, r093, r094, r095, r096, r097)
    b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116, r117)
    b12Θ = @evalpoly(Θ, 0, 0, r122, r123, r124, r125, r126, r127)
    b13Θ = @evalpoly(Θ, 0, 0, r132, r133, r134, r135, r136, r137)
    b14Θ = @evalpoly(Θ, 0, 0, r142, r143, r144, r145, r146, r147)
    b15Θ = @evalpoly(Θ, 0, 0, r152, r153, r154, r155, r156, r157)
    b16Θ = @evalpoly(Θ, 0, 0, r162, r163, r164, r165, r166, r167)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern7ConstantCache,
                                  idxs::Nothing, T::Type{Val{0}})
    @vern7pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[11]*b11Θ + k[12]*b12Θ + k[13]*b13Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ + k[7] * b7Θ +
                      k[8] * b8Θ + k[9] * b9Θ + k[11] * b11Θ + k[12] * b12Θ + k[13] * b13Θ +
                      k[14] * b14Θ + k[15] * b15Θ + k[16] * b16Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern7Cache, idxs::Nothing,
                                  T::Type{Val{0}})
    @vern7pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[4]*b4Θ + k[5]*b5Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[11]*b11Θ + k[12]*b12Θ + k[13]*b13Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ)
    return @inbounds @.. broadcast=false y₀+dt * (k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ +
                                             k[6] * b6Θ + k[7] * b7Θ +
                                             k[8] * b8Θ + k[9] * b9Θ + k[11] * b11Θ +
                                             k[12] * b12Θ + k[13] * b13Θ +
                                             k[14] * b14Θ + k[15] * b15Θ + k[16] * b16Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern7ConstantCache, Vern7Cache}, idxs,
                                  T::Type{Val{0}})
    @vern7pre0
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
            k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
            k[11][idxs] * b11Θ + k[12][idxs] * b12Θ + k[13][idxs] * b13Θ +
            k[14][idxs] * b14Θ + k[15][idxs] * b15Θ + k[16][idxs] * b16Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern7ConstantCache, Vern7Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @vern7pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[4] * b4Θ + k[5] * b5Θ + k[6] * b6Θ +
                                       k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ + k[11] * b11Θ +
                                       k[12] * b12Θ + k[13] * b13Θ + k[14] * b14Θ +
                                       k[15] * b15Θ + k[16] * b16Θ)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern7ConstantCache, Vern7Cache}, idxs,
                                   T::Type{Val{0}})
    @vern7pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[4][idxs] * b4Θ + k[5][idxs] * b5Θ +
                                    k[6][idxs] * b6Θ + k[7][idxs] * b7Θ + k[8][idxs] * b8Θ +
                                    k[9][idxs] * b9Θ + k[11][idxs] * b11Θ +
                                    k[12][idxs] * b12Θ + k[13][idxs] * b13Θ +
                                    k[14][idxs] * b14Θ + k[15][idxs] * b15Θ +
                                    k[16][idxs] * b16Θ)
    out
end

@def vern7pre1 begin
    @vern7unpack
    b1Θdiff = @evalpoly(Θ, r011, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016, 7*r017)
    b4Θdiff = @evalpoly(Θ, 0, 2*r042, 3*r043, 4*r044, 5*r045, 6*r046, 7*r047)
    b5Θdiff = @evalpoly(Θ, 0, 2*r052, 3*r053, 4*r054, 5*r055, 6*r056, 7*r057)
    b6Θdiff = @evalpoly(Θ, 0, 2*r062, 3*r063, 4*r064, 5*r065, 6*r066, 7*r067)
    b7Θdiff = @evalpoly(Θ, 0, 2*r072, 3*r073, 4*r074, 5*r075, 6*r076, 7*r077)
    b8Θdiff = @evalpoly(Θ, 0, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086, 7*r087)
    b9Θdiff = @evalpoly(Θ, 0, 2*r092, 3*r093, 4*r094, 5*r095, 6*r096, 7*r097)
    b11Θdiff = @evalpoly(Θ, 0, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116, 7*r117)
    b12Θdiff = @evalpoly(Θ, 0, 2*r122, 3*r123, 4*r124, 5*r125, 6*r126, 7*r127)
    b13Θdiff = @evalpoly(Θ, 0, 2*r132, 3*r133, 4*r134, 5*r135, 6*r136, 7*r137)
    b14Θdiff = @evalpoly(Θ, 0, 2*r142, 3*r143, 4*r144, 5*r145, 6*r146, 7*r147)
    b15Θdiff = @evalpoly(Θ, 0, 2*r152, 3*r153, 4*r154, 5*r155, 6*r156, 7*r157)
    b16Θdiff = @evalpoly(Θ, 0, 2*r162, 3*r163, 4*r164, 5*r165, 6*r166, 7*r167)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern7ConstantCache, Vern7Cache},
                                  idxs::Nothing, T::Type{Val{1}})
    @vern7pre1
    #@.. broadcast=false k[1]*b1Θdiff + k[4]*b4Θdiff + k[5]*b5Θdiff + k[6]*b6Θdiff + k[7]*b7Θdiff + k[8]*b8Θdiff + k[9]*b9Θdiff + k[11]*b11Θdiff + k[12]*b12Θdiff + k[13]*b13Θdiff + k[14]*b14Θdiff + k[15]*b15Θdiff + k[16]*b16Θdiff
    return @inbounds k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff + k[6] * b6Θdiff +
                     k[7] * b7Θdiff + k[8] * b8Θdiff + k[9] * b9Θdiff + k[11] * b11Θdiff +
                     k[12] * b12Θdiff + k[13] * b13Θdiff + k[14] * b14Θdiff +
                     k[15] * b15Θdiff +
                     k[16] * b16Θdiff
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern7ConstantCache, Vern7Cache}, idxs,
                                  T::Type{Val{1}})
    @vern7pre1
    return k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff + k[5][idxs] * b5Θdiff +
           k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff + k[8][idxs] * b8Θdiff +
           k[9][idxs] * b9Θdiff + k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff +
           k[13][idxs] * b13Θdiff + k[14][idxs] * b14Θdiff + k[15][idxs] * b15Θdiff +
           k[16][idxs] * b16Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern7ConstantCache, Vern7Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @vern7pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[4] * b4Θdiff + k[5] * b5Θdiff +
                                      k[6] * b6Θdiff + k[7] * b7Θdiff + k[8] * b8Θdiff +
                                      k[9] * b9Θdiff + k[11] * b11Θdiff + k[12] * b12Θdiff +
                                      k[13] * b13Θdiff + k[14] * b14Θdiff +
                                      k[15] * b15Θdiff + k[16] * b16Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern7ConstantCache, Vern7Cache}, idxs,
                                   T::Type{Val{1}})
    @vern7pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[4][idxs] * b4Θdiff +
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
    @unpack r011, r012, r013, r014, r015, r016, r017, r018, r062, r063, r064, r065, r066, r067, r068, r072, r073, r074, r075, r076, r077, r078, r082, r083, r084, r085, r086, r087, r088, r092, r093, r094, r095, r096, r097, r098, r102, r103, r104, r105, r106, r107, r108, r112, r113, r114, r115, r116, r117, r118, r122, r123, r124, r125, r126, r127, r128, r142, r143, r144, r145, r146, r147, r148, r152, r153, r154, r155, r156, r157, r158, r162, r163, r164, r165, r166, r167, r168, r172, r173, r174, r175, r176, r177, r178, r182, r183, r184, r185, r186, r187, r188, r192, r193, r194, r195, r196, r197, r198, r202, r203, r204, r205, r206, r207, r208, r212, r213, r214, r215, r216, r217, r218 = cache.tab.interp
end

@def vern8pre0 begin
    @vern8unpack
    b1Θ = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017, r018)
    b6Θ = @evalpoly(Θ, 0, 0, r062, r063, r064, r065, r066, r067, r068)
    b7Θ = @evalpoly(Θ, 0, 0, r072, r073, r074, r075, r076, r077, r078)
    b8Θ = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086, r087, r088)
    b9Θ = @evalpoly(Θ, 0, 0, r092, r093, r094, r095, r096, r097, r098)
    b10Θ = @evalpoly(Θ, 0, 0, r102, r103, r104, r105, r106, r107, r108)
    b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116, r117, r118)
    b12Θ = @evalpoly(Θ, 0, 0, r122, r123, r124, r125, r126, r127, r128)
    b14Θ = @evalpoly(Θ, 0, 0, r142, r143, r144, r145, r146, r147, r148)
    b15Θ = @evalpoly(Θ, 0, 0, r152, r153, r154, r155, r156, r157, r158)
    b16Θ = @evalpoly(Θ, 0, 0, r162, r163, r164, r165, r166, r167, r168)
    b17Θ = @evalpoly(Θ, 0, 0, r172, r173, r174, r175, r176, r177, r178)
    b18Θ = @evalpoly(Θ, 0, 0, r182, r183, r184, r185, r186, r187, r188)
    b19Θ = @evalpoly(Θ, 0, 0, r192, r193, r194, r195, r196, r197, r198)
    b20Θ = @evalpoly(Θ, 0, 0, r202, r203, r204, r205, r206, r207, r208)
    b21Θ = @evalpoly(Θ, 0, 0, r212, r213, r214, r215, r216, r217, r218)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern8ConstantCache,
                                  idxs::Nothing, T::Type{Val{0}})
    @vern8pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ + k[17]*b17Θ + k[18]*b18Θ + k[19]*b19Θ + k[20]*b20Θ + k[21]*b21Θ)
    return @inbounds y₀ +
                     dt * (k[1] * b1Θ + k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ + k[9] * b9Θ +
                      k[10] * b10Θ + k[11] * b11Θ + k[12] * b12Θ + k[14] * b14Θ +
                      k[15] * b15Θ +
                      k[16] * b16Θ + k[17] * b17Θ + k[18] * b18Θ + k[19] * b19Θ +
                      k[20] * b20Θ +
                      k[21] * b21Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern8Cache, idxs::Nothing,
                                  T::Type{Val{0}})
    @vern8pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[6]*b6Θ + k[7]*b7Θ + k[8]*b8Θ + k[9]*b9Θ + k[10]*b10Θ + k[11]*b11Θ + k[12]*b12Θ + k[14]*b14Θ + k[15]*b15Θ + k[16]*b16Θ + k[17]*b17Θ + k[18]*b18Θ + k[19]*b19Θ + k[20]*b20Θ + k[21]*b21Θ)
    return @inbounds @.. broadcast=false y₀+dt * (k[1] * b1Θ + k[6] * b6Θ + k[7] * b7Θ +
                                             k[8] * b8Θ + k[9] * b9Θ +
                                             k[10] * b10Θ + k[11] * b11Θ + k[12] * b12Θ +
                                             k[14] * b14Θ + k[15] * b15Θ +
                                             k[16] * b16Θ + k[17] * b17Θ + k[18] * b18Θ +
                                             k[19] * b19Θ + k[20] * b20Θ +
                                             k[21] * b21Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
                                  T::Type{Val{0}})
    @vern8pre0
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ +
            k[8][idxs] * b8Θ + k[9][idxs] * b9Θ + k[10][idxs] * b10Θ +
            k[11][idxs] * b11Θ + k[12][idxs] * b12Θ + k[14][idxs] * b14Θ +
            k[15][idxs] * b15Θ + k[16][idxs] * b16Θ + k[17][idxs] * b17Θ +
            k[18][idxs] * b18Θ + k[19][idxs] * b19Θ + k[20][idxs] * b20Θ +
            k[21][idxs] * b21Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern8ConstantCache, Vern8Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @vern8pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[6] * b6Θ + k[7] * b7Θ + k[8] * b8Θ +
                                       k[9] * b9Θ + k[10] * b10Θ + k[11] * b11Θ +
                                       k[12] * b12Θ + k[14] * b14Θ + k[15] * b15Θ +
                                       k[16] * b16Θ + k[17] * b17Θ + k[18] * b18Θ +
                                       k[19] * b19Θ + k[20] * b20Θ + k[21] * b21Θ)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
                                   T::Type{Val{0}})
    @vern8pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[6][idxs] * b6Θ + k[7][idxs] * b7Θ +
                                    k[8][idxs] * b8Θ + k[9][idxs] * b9Θ +
                                    k[10][idxs] * b10Θ + k[11][idxs] * b11Θ +
                                    k[12][idxs] * b12Θ + k[14][idxs] * b14Θ +
                                    k[15][idxs] * b15Θ + k[16][idxs] * b16Θ +
                                    k[17][idxs] * b17Θ + k[18][idxs] * b18Θ +
                                    k[19][idxs] * b19Θ + k[20][idxs] * b20Θ +
                                    k[21][idxs] * b21Θ)
    out
end

@def vern8pre1 begin
    @vern8unpack
    b1Θdiff = @evalpoly(Θ, r011, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016, 7*r017, 8*r018)
    b6Θdiff = @evalpoly(Θ, 0, 2*r062, 3*r063, 4*r064, 5*r065, 6*r066, 7*r067, 8*r068)
    b7Θdiff = @evalpoly(Θ, 0, 2*r072, 3*r073, 4*r074, 5*r075, 6*r076, 7*r077, 8*r078)
    b8Θdiff = @evalpoly(Θ, 0, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086, 7*r087, 8*r088)
    b9Θdiff = @evalpoly(Θ, 0, 2*r092, 3*r093, 4*r094, 5*r095, 6*r096, 7*r097, 8*r098)
    b10Θdiff = @evalpoly(Θ, 0, 2*r102, 3*r103, 4*r104, 5*r105, 6*r106, 7*r107, 8*r108)
    b11Θdiff = @evalpoly(Θ, 0, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116, 7*r117, 8*r118)
    b12Θdiff = @evalpoly(Θ, 0, 2*r122, 3*r123, 4*r124, 5*r125, 6*r126, 7*r127, 8*r128)
    b14Θdiff = @evalpoly(Θ, 0, 2*r142, 3*r143, 4*r144, 5*r145, 6*r146, 7*r147, 8*r148)
    b15Θdiff = @evalpoly(Θ, 0, 2*r152, 3*r153, 4*r154, 5*r155, 6*r156, 7*r157, 8*r158)
    b16Θdiff = @evalpoly(Θ, 0, 2*r162, 3*r163, 4*r164, 5*r165, 6*r166, 7*r167, 8*r168)
    b17Θdiff = @evalpoly(Θ, 0, 2*r172, 3*r173, 4*r174, 5*r175, 6*r176, 7*r177, 8*r178)
    b18Θdiff = @evalpoly(Θ, 0, 2*r182, 3*r183, 4*r184, 5*r185, 6*r186, 7*r187, 8*r188)
    b19Θdiff = @evalpoly(Θ, 0, 2*r192, 3*r193, 4*r194, 5*r195, 6*r196, 7*r197, 8*r198)
    b20Θdiff = @evalpoly(Θ, 0, 2*r202, 3*r203, 4*r204, 5*r205, 6*r206, 7*r207, 8*r208)
    b21Θdiff = @evalpoly(Θ, 0, 2*r212, 3*r213, 4*r214, 5*r215, 6*r216, 7*r217, 8*r218)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern8ConstantCache, Vern8Cache},
                                  idxs::Nothing, T::Type{Val{1}})
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

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
                                  T::Type{Val{1}})
    @vern8pre1
    return k[1][idxs] * b1Θdiff + k[6][idxs] * b6Θdiff + k[7][idxs] * b7Θdiff +
           k[8][idxs] * b8Θdiff + k[9][idxs] * b9Θdiff + k[10][idxs] * b10Θdiff +
           k[11][idxs] * b11Θdiff + k[12][idxs] * b12Θdiff + k[14][idxs] * b14Θdiff +
           k[15][idxs] * b15Θdiff + k[16][idxs] * b16Θdiff + k[17][idxs] * b17Θdiff +
           k[18][idxs] * b18Θdiff + k[19][idxs] * b19Θdiff + k[20][idxs] * b20Θdiff +
           k[21][idxs] * b21Θdiff
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern8ConstantCache, Vern8Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @vern8pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[6] * b6Θdiff + k[7] * b7Θdiff +
                                      k[8] * b8Θdiff + k[9] * b9Θdiff + k[10] * b10Θdiff +
                                      k[11] * b11Θdiff + k[12] * b12Θdiff +
                                      k[14] * b14Θdiff + k[15] * b15Θdiff +
                                      k[16] * b16Θdiff + k[17] * b17Θdiff +
                                      k[18] * b18Θdiff + k[19] * b19Θdiff +
                                      k[20] * b20Θdiff + k[21] * b21Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern8ConstantCache, Vern8Cache}, idxs,
                                   T::Type{Val{1}})
    @vern8pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[6][idxs] * b6Θdiff +
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
    T = Float64
    r011 = convert(T, 1)
    r012 = convert(T, -28.330488700617398)
    r013 = convert(T, 257.6535452078578)
    r014 = convert(T, -1152.1544557434572)
    r015 = convert(T, 2909.390878345409)
    r016 = convert(T, -4355.005172868188)
    r017 = convert(T, 3834.083497036262)
    r018 = convert(T, -1835.419052683407)
    r019 = convert(T, 368.7958613829998)
    r082 = convert(T, 2.649656243770091)
    r083 = convert(T, -96.30312807816006)
    r084 = convert(T, 869.3095462492796)
    r085 = convert(T, -3395.688567551074)
    r086 = convert(T, 6796.933987158715)
    r087 = convert(T, -7340.848417712072)
    r088 = convert(T, 4082.8488969923656)
    r089 = convert(T, -919.2934944890586)
    r092 = convert(T, -1.5639451819287329)
    r093 = convert(T, 56.8423973927286)
    r094 = convert(T, -513.1052300304285)
    r095 = convert(T, 2004.2867021103232)
    r096 = convert(T, -4011.8533059139295)
    r097 = convert(T, 4332.895839278586)
    r098 = convert(T, -2409.8793479371448)
    r099 = convert(T, 542.6079835318221)
    r102 = convert(T, -0.8627103334967224)
    r103 = convert(T, 31.355653751851733)
    r104 = convert(T, -283.0413682227354)
    r105 = convert(T, 1105.613463426007)
    r106 = convert(T, -2213.0362006784526)
    r107 = convert(T, 2390.1310977541207)
    r108 = convert(T, -1329.3482661468738)
    r109 = convert(T, 299.31580712657853)
    r112 = convert(T, -1.5202953379012147)
    r113 = convert(T, 55.25592121120227)
    r114 = convert(T, -498.7844190970741)
    r115 = convert(T, 1948.346888525776)
    r116 = convert(T, -3899.8821364075516)
    r117 = convert(T, 4211.964345158858)
    r118 = convert(T, -2342.619408856117)
    r119 = convert(T, 527.4637482204279)
    r122 = convert(T, -3.8469388441255234)
    r123 = convert(T, 139.81898409868404)
    r124 = convert(T, -1262.1186876216004)
    r125 = convert(T, 4930.075848057311)
    r126 = convert(T, -9868.21948606954)
    r127 = convert(T, 10657.908924348867)
    r128 = convert(T, -5927.738759872814)
    r129 = convert(T, 1334.688551172191)
    r132 = convert(T, -0.39427130612001415)
    r133 = convert(T, 14.329994760676497)
    r134 = convert(T, -129.35406659945582)
    r135 = convert(T, 505.28160770025175)
    r136 = convert(T, -1011.3900801394333)
    r137 = convert(T, 1092.3250517818917)
    r138 = convert(T, -607.531701930281)
    r139 = convert(T, 136.79172444804232)
    r142 = convert(T, -0.9233145622082102)
    r143 = convert(T, 33.55834582309799)
    r144 = convert(T, -302.9246397549736)
    r145 = convert(T, 1183.2813069678675)
    r146 = convert(T, -2368.4989867901113)
    r147 = convert(T, 2558.034559755808)
    r148 = convert(T, -1422.7331755778803)
    r149 = convert(T, 320.3423358787482)
    r152 = convert(T, -0.20688628029300538)
    r153 = convert(T, 7.519388975651663)
    r154 = convert(T, -67.87605708082904)
    r155 = convert(T, 265.136799698415)
    r156 = convert(T, -530.7074807559026)
    r157 = convert(T, 573.176549564149)
    r158 = convert(T, -318.7905688834869)
    r159 = convert(T, 71.77882490212657)
    r172 = convert(T, -0.44724419067440996)
    r173 = convert(T, 16.44684676010504)
    r174 = convert(T, -154.40861059212955)
    r175 = convert(T, 641.8986298540249)
    r176 = convert(T, -1391.9392256879823)
    r177 = convert(T, 1643.890568302952)
    r178 = convert(T, -1004.0652972233179)
    r179 = convert(T, 248.6243327770223)
    r182 = convert(T, -0.1507876007899798)
    r183 = convert(T, 5.527328824824632)
    r184 = convert(T, -51.33833743084619)
    r185 = convert(T, 209.60220027032804)
    r186 = convert(T, -442.7692650421826)
    r187 = convert(T, 505.0579312588053)
    r188 = convert(T, -295.63364106156195)
    r189 = convert(T, 69.70457078142275)
    r192 = convert(T, -0.6413652207435296)
    r193 = convert(T, 23.510132486246846)
    r194 = convert(T, -218.36426832469724)
    r195 = convert(T, 891.5292818535365)
    r196 = convert(T, -1883.290177206008)
    r197 = convert(T, 2148.2309544883997)
    r198 = convert(T, -1257.4584015217124)
    r199 = convert(T, 296.4838434449778)
    r202 = convert(T, 1.8107293134448457)
    r203 = convert(T, -66.37479657295337)
    r204 = convert(T, 616.4952025401107)
    r205 = convert(T, -2517.0030307773227)
    r206 = convert(T, 5316.984175781034)
    r207 = convert(T, -6064.976140789574)
    r208 = convert(T, 3550.1095388883914)
    r209 = convert(T, -837.0456783831302)
    r212 = convert(T, 0.05176008760353718)
    r213 = convert(T, -1.8973378625803488)
    r214 = convert(T, 17.622648207936294)
    r215 = convert(T, -71.94907400242467)
    r216 = convert(T, 151.9871383765666)
    r217 = convert(T, -173.36864987478606)
    r218 = convert(T, 101.4806461521468)
    r219 = convert(T, -23.927131084462175)
    r222 = convert(T, 31.321782556688)
    r223 = convert(T, -355.6570858339106)
    r224 = convert(T, 1752.6852824895159)
    r225 = convert(T, -4708.092293138363)
    r226 = convert(T, 7370.900776193489)
    r227 = convert(T, -6716.504964764566)
    r228 = convert(T, 3303.940398161186)
    r229 = convert(T, -678.5938956640391)
    r232 = convert(T, -2.7196073341859246)
    r233 = convert(T, 86.64045615858264)
    r234 = convert(T, -454.1926030939031)
    r235 = convert(T, 1014.7492211005434)
    r236 = convert(T, -1133.583456714544)
    r237 = convert(T, 610.4671827718666)
    r238 = convert(T, -109.02334994495438)
    r239 = convert(T, -12.337842943405471)
    r242 = convert(T, 3.1772148014329233)
    r243 = convert(T, -113.8098697715143)
    r244 = convert(T, 978.0935981825675)
    r245 = convert(T, -3575.1293776236703)
    r246 = convert(T, 6764.3615198384505)
    r247 = convert(T, -6987.161043852012)
    r248 = convert(T, 3751.9057627895713)
    r249 = convert(T, -821.4378043648254)
    r252 = convert(T, 0.877284308346553)
    r253 = convert(T, -31.51810423988375)
    r254 = convert(T, 273.1229151353221)
    r255 = convert(T, -993.2198643101782)
    r256 = convert(T, 1787.888078312664)
    r257 = convert(T, -1677.394835799641)
    r258 = convert(T, 781.3579535062688)
    r259 = convert(T, -141.11342691289855)
    r262 = convert(T, 1.7194275817987157)
    r263 = convert(T, -62.89867309250732)
    r264 = convert(T, 580.333550787398)
    r265 = convert(T, -2348.110620506761)
    r266 = convert(T, 4921.119298612906)
    r267 = convert(T, -5597.912448707917)
    r268 = convert(T, 3288.5977751496216)
    r269 = convert(T, -782.8483098245397)
end

@def vern9pre0 begin
    @vern9unpack
    b1Θ = @evalpoly(Θ, 0, r011, r012, r013, r014, r015, r016, r017, r018, r019)
    b8Θ = @evalpoly(Θ, 0, 0, r082, r083, r084, r085, r086, r087, r088, r089)
    b9Θ = @evalpoly(Θ, 0, 0, r092, r093, r094, r095, r096, r097, r098, r099)
    b10Θ = @evalpoly(Θ, 0, 0, r102, r103, r104, r105, r106, r107, r108, r109)
    b11Θ = @evalpoly(Θ, 0, 0, r112, r113, r114, r115, r116, r117, r118, r119)
    b12Θ = @evalpoly(Θ, 0, 0, r122, r123, r124, r125, r126, r127, r128, r129)
    b13Θ = @evalpoly(Θ, 0, 0, r132, r133, r134, r135, r136, r137, r138, r139)
    b14Θ = @evalpoly(Θ, 0, 0, r142, r143, r144, r145, r146, r147, r148, r149)
    b15Θ = @evalpoly(Θ, 0, 0, r152, r153, r154, r155, r156, r157, r158, r159)
    b17Θ = @evalpoly(Θ, 0, 0, r172, r173, r174, r175, r176, r177, r178, r179)
    b18Θ = @evalpoly(Θ, 0, 0, r182, r183, r184, r185, r186, r187, r188, r189)
    b19Θ = @evalpoly(Θ, 0, 0, r192, r193, r194, r195, r196, r197, r198, r199)
    b20Θ = @evalpoly(Θ, 0, 0, r202, r203, r204, r205, r206, r207, r208, r209)
    b21Θ = @evalpoly(Θ, 0, 0, r212, r213, r214, r215, r216, r217, r218, r219)
    b22Θ = @evalpoly(Θ, 0, 0, r222, r223, r224, r225, r226, r227, r228, r229)
    b23Θ = @evalpoly(Θ, 0, 0, r232, r233, r234, r235, r236, r237, r238, r239)
    b24Θ = @evalpoly(Θ, 0, 0, r242, r243, r244, r245, r246, r247, r248, r249)
    b25Θ = @evalpoly(Θ, 0, 0, r252, r253, r254, r255, r256, r257, r258, r259)
    b26Θ = @evalpoly(Θ, 0, 0, r262, r263, r264, r265, r266, r267, r268, r269)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern9ConstantCache,
                                  idxs::Nothing, T::Type{Val{0}})
    @vern9pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[2]*b8Θ + k[3]*b9Θ + k[4]*b10Θ + k[5]*b11Θ + k[6]*b12Θ + k[7]*b13Θ + k[8]*b14Θ + k[9]*b15Θ + k[11]*b17Θ + k[12]*b18Θ + k[13]*b19Θ + k[14]*b20Θ + k[15]*b21Θ + k[16]*b22Θ + k[17]*b23Θ + k[18]*b24Θ + k[19]*b25Θ + k[20]*b26Θ)
    return @inbounds y₀ +
                     dt *
                     (k[1] * b1Θ + k[2] * b8Θ + k[3] * b9Θ + k[4] * b10Θ + k[5] * b11Θ +
                      k[6] * b12Θ + k[7] * b13Θ + k[8] * b14Θ + k[9] * b15Θ + k[11] * b17Θ +
                      k[12] * b18Θ + k[13] * b19Θ + k[14] * b20Θ + k[15] * b21Θ +
                      k[16] * b22Θ +
                      k[17] * b23Θ + k[18] * b24Θ + k[19] * b25Θ + k[20] * b26Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k, cache::Vern9Cache, idxs::Nothing,
                                  T::Type{Val{0}})
    @vern9pre0
    #@.. broadcast=false y₀ + dt*(k[1]*b1Θ + k[2]*b8Θ + k[3]*b9Θ + k[4]*b10Θ + k[5]*b11Θ + k[6]*b12Θ + k[7]*b13Θ + k[8]*b14Θ + k[9]*b15Θ + k[11]*b17Θ + k[12]*b18Θ + k[13]*b19Θ + k[14]*b20Θ + k[15]*b21Θ + k[16]*b22Θ + k[17]*b23Θ + k[18]*b24Θ + k[19]*b25Θ + k[20]*b26Θ)
    return @inbounds @.. broadcast=false y₀+dt * (k[1] * b1Θ + k[2] * b8Θ + k[3] * b9Θ +
                                             k[4] * b10Θ + k[5] * b11Θ +
                                             k[6] * b12Θ + k[7] * b13Θ + k[8] * b14Θ +
                                             k[9] * b15Θ + k[11] * b17Θ +
                                             k[12] * b18Θ + k[13] * b19Θ + k[14] * b20Θ +
                                             k[15] * b21Θ + k[16] * b22Θ +
                                             k[17] * b23Θ + k[18] * b24Θ + k[19] * b25Θ +
                                             k[20] * b26Θ)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                  T::Type{Val{0}})
    @vern9pre0
    return y₀[idxs] +
           dt * (k[1][idxs] * b1Θ + k[2][idxs] * b8Θ + k[3][idxs] * b9Θ +
            k[4][idxs] * b10Θ + k[5][idxs] * b11Θ + k[6][idxs] * b12Θ +
            k[7][idxs] * b13Θ + k[8][idxs] * b14Θ + k[9][idxs] * b15Θ + k[11][idxs] * b17Θ +
            k[12][idxs] * b18Θ + k[13][idxs] * b19Θ + k[14][idxs] * b20Θ +
            k[15][idxs] * b21Θ +
            k[16][idxs] * b22Θ + k[17][idxs] * b23Θ + k[18][idxs] * b24Θ +
            k[19][idxs] * b25Θ + k[20][idxs] * b26Θ)
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @vern9pre0
    @inbounds @.. broadcast=false out=y₀ +
                                      dt *
                                      (k[1] * b1Θ + k[2] * b8Θ + k[3] * b9Θ + k[4] * b10Θ +
                                       k[5] * b11Θ + k[6] * b12Θ + k[7] * b13Θ +
                                       k[8] * b14Θ + k[9] * b15Θ + k[11] * b17Θ +
                                       k[12] * b18Θ + k[13] * b19Θ + k[14] * b20Θ +
                                       k[15] * b21Θ + k[16] * b22Θ + k[17] * b23Θ +
                                       k[18] * b24Θ + k[19] * b25Θ + k[20] * b26Θ)
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                   T::Type{Val{0}})
    @vern9pre0
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt *
                                   (k[1][idxs] * b1Θ + k[2][idxs] * b8Θ + k[3][idxs] * b9Θ +
                                    k[4][idxs] * b10Θ + k[5][idxs] * b11Θ +
                                    k[6][idxs] * b12Θ + k[7][idxs] * b13Θ +
                                    k[8][idxs] * b14Θ + k[9][idxs] * b15Θ +
                                    k[11][idxs] * b17Θ + k[12][idxs] * b18Θ +
                                    k[13][idxs] * b19Θ + k[14][idxs] * b20Θ +
                                    k[15][idxs] * b21Θ + k[16][idxs] * b22Θ +
                                    k[17][idxs] * b23Θ + k[18][idxs] * b24Θ +
                                    k[19][idxs] * b25Θ + k[20][idxs] * b26Θ)
    out
end

@def vern9pre1 begin
    @vern9unpack
    b1Θdiff = @evalpoly(Θ, r011, 2*r012, 3*r013, 4*r014, 5*r015, 6*r016, 7*r017, 8*r018,
                        9*r019)
    b8Θdiff = @evalpoly(Θ, 0, 2*r082, 3*r083, 4*r084, 5*r085, 6*r086, 7*r087, 8*r088,
                        9*r089)
    b9Θdiff = @evalpoly(Θ, 0, 2*r092, 3*r093, 4*r094, 5*r095, 6*r096, 7*r097, 8*r098,
                        9*r099)
    b10Θdiff = @evalpoly(Θ, 0, 2*r102, 3*r103, 4*r104, 5*r105, 6*r106, 7*r107, 8*r108,
                         9*r109)
    b11Θdiff = @evalpoly(Θ, 0, 2*r112, 3*r113, 4*r114, 5*r115, 6*r116, 7*r117, 8*r118,
                         9*r119)
    b12Θdiff = @evalpoly(Θ, 0, 2*r122, 3*r123, 4*r124, 5*r125, 6*r126, 7*r127, 8*r128,
                         9*r129)
    b13Θdiff = @evalpoly(Θ, 0, 2*r132, 3*r133, 4*r134, 5*r135, 6*r136, 7*r137, 8*r138,
                         9*r139)
    b14Θdiff = @evalpoly(Θ, 0, 2*r142, 3*r143, 4*r144, 5*r145, 6*r146, 7*r147, 8*r148,
                         9*r149)
    b15Θdiff = @evalpoly(Θ, 0, 2*r152, 3*r153, 4*r154, 5*r155, 6*r156, 7*r157, 8*r158,
                         9*r159)
    b17Θdiff = @evalpoly(Θ, 0, 2*r172, 3*r173, 4*r174, 5*r175, 6*r176, 7*r177, 8*r178,
                         9*r179)
    b18Θdiff = @evalpoly(Θ, 0, 2*r182, 3*r183, 4*r184, 5*r185, 6*r186, 7*r187, 8*r188,
                         9*r189)
    b19Θdiff = @evalpoly(Θ, 0, 2*r192, 3*r193, 4*r194, 5*r195, 6*r196, 7*r197, 8*r198,
                         9*r199)
    b20Θdiff = @evalpoly(Θ, 0, 2*r202, 3*r203, 4*r204, 5*r205, 6*r206, 7*r207, 8*r208,
                         9*r209)
    b21Θdiff = @evalpoly(Θ, 0, 2*r212, 3*r213, 4*r214, 5*r215, 6*r216, 7*r217, 8*r218,
                         9*r219)
    b22Θdiff = @evalpoly(Θ, 0, 2*r222, 3*r223, 4*r224, 5*r225, 6*r226, 7*r227, 8*r228,
                         9*r229)
    b23Θdiff = @evalpoly(Θ, 0, 2*r232, 3*r233, 4*r234, 5*r235, 6*r236, 7*r237, 8*r238,
                         9*r239)
    b24Θdiff = @evalpoly(Θ, 0, 2*r242, 3*r243, 4*r244, 5*r245, 6*r246, 7*r247, 8*r248,
                         9*r249)
    b25Θdiff = @evalpoly(Θ, 0, 2*r252, 3*r253, 4*r254, 5*r255, 6*r256, 7*r257, 8*r258,
                         9*r259)
    b26Θdiff = @evalpoly(Θ, 0, 2*r262, 3*r263, 4*r264, 5*r265, 6*r266, 7*r267, 8*r268,
                         9*r269)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache},
                                  idxs::Nothing, T::Type{Val{1}})
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

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                  T::Type{Val{1}})
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

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache},
                                   idxs::Nothing, T::Type{Val{1}})
    @vern9pre1
    @inbounds @.. broadcast=false out=k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
                                      k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
                                      k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
                                      k[11] * b17Θdiff + k[12] * b18Θdiff +
                                      k[13] * b19Θdiff + k[14] * b20Θdiff +
                                      k[15] * b21Θdiff + k[16] * b22Θdiff +
                                      k[17] * b23Θdiff + k[18] * b24Θdiff +
                                      k[19] * b25Θdiff + k[20] * b26Θdiff
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                   T::Type{Val{1}})
    @vern9pre1
    @views @.. broadcast=false out=k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
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
    b1Θdiff = @evalpoly(Θ, 2*r012, 6*r013, 12*r014, 20*r015, 30*r016, 42*r017, 56*r018,
                        72*r019)
    b8Θdiff = @evalpoly(Θ, 2*r082, 6*r083, 12*r084, 20*r085, 30*r086, 42*r087, 56*r088,
                        72*r089)
    b9Θdiff = @evalpoly(Θ, 2*r092, 6*r093, 12*r094, 20*r095, 30*r096, 42*r097, 56*r098,
                        72*r099)
    b10Θdiff = @evalpoly(Θ, 2*r102, 6*r103, 12*r104, 20*r105, 30*r106, 42*r107, 56*r108,
                         72*r109)
    b11Θdiff = @evalpoly(Θ, 2*r112, 6*r113, 12*r114, 20*r115, 30*r116, 42*r117, 56*r118,
                         72*r119)
    b12Θdiff = @evalpoly(Θ, 2*r122, 6*r123, 12*r124, 20*r125, 30*r126, 42*r127, 56*r128,
                         72*r129)
    b13Θdiff = @evalpoly(Θ, 2*r132, 6*r133, 12*r134, 20*r135, 30*r136, 42*r137, 56*r138,
                         72*r139)
    b14Θdiff = @evalpoly(Θ, 2*r142, 6*r143, 12*r144, 20*r145, 30*r146, 42*r147, 56*r148,
                         72*r149)
    b15Θdiff = @evalpoly(Θ, 2*r152, 6*r153, 12*r154, 20*r155, 30*r156, 42*r157, 56*r158,
                         72*r159)
    b17Θdiff = @evalpoly(Θ, 2*r172, 6*r173, 12*r174, 20*r175, 30*r176, 42*r177, 56*r178,
                         72*r179)
    b18Θdiff = @evalpoly(Θ, 2*r182, 6*r183, 12*r184, 20*r185, 30*r186, 42*r187, 56*r188,
                         72*r189)
    b19Θdiff = @evalpoly(Θ, 2*r192, 6*r193, 12*r194, 20*r195, 30*r196, 42*r197, 56*r198,
                         72*r199)
    b20Θdiff = @evalpoly(Θ, 2*r202, 6*r203, 12*r204, 20*r205, 30*r206, 42*r207, 56*r208,
                         72*r209)
    b21Θdiff = @evalpoly(Θ, 2*r212, 6*r213, 12*r214, 20*r215, 30*r216, 42*r217, 56*r218,
                         72*r219)
    b22Θdiff = @evalpoly(Θ, 2*r222, 6*r223, 12*r224, 20*r225, 30*r226, 42*r227, 56*r228,
                         72*r229)
    b23Θdiff = @evalpoly(Θ, 2*r232, 6*r233, 12*r234, 20*r235, 30*r236, 42*r237, 56*r238,
                         72*r239)
    b24Θdiff = @evalpoly(Θ, 2*r242, 6*r243, 12*r244, 20*r245, 30*r246, 42*r247, 56*r248,
                         72*r249)
    b25Θdiff = @evalpoly(Θ, 2*r252, 6*r253, 12*r254, 20*r255, 30*r256, 42*r257, 56*r258,
                         72*r259)
    b26Θdiff = @evalpoly(Θ, 2*r262, 6*r263, 12*r264, 20*r265, 30*r266, 42*r267, 56*r268,
                         72*r269)
    invdt = inv(dt)
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache},
                                  idxs::Nothing, T::Type{Val{2}})
    @vern9pre2
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds (k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
                      k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff +
                      k[8] * b14Θdiff +
                      k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
                      k[13] * b19Θdiff +
                      k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
                      k[17] * b23Θdiff +
                      k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff) * invdt
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                  T::Type{Val{2}})
    @vern9pre2
    return (k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
            k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
            k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
            k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
            k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
            k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
            k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
            k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
            k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff) * invdt
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache},
                                   idxs::Nothing, T::Type{Val{2}})
    @vern9pre2
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
                                       k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
                                       k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
                                       k[11] * b17Θdiff + k[12] * b18Θdiff +
                                       k[13] * b19Θdiff + k[14] * b20Θdiff +
                                       k[15] * b21Θdiff + k[16] * b22Θdiff +
                                       k[17] * b23Θdiff + k[18] * b24Θdiff +
                                       k[19] * b25Θdiff + k[20] * b26Θdiff) * invdt
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                   T::Type{Val{2}})
    @vern9pre2
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
                                    k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
                                    k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
                                    k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
                                    k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
                                    k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
                                    k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
                                    k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
                                    k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
                                    k[20][idxs] * b26Θdiff) * invdt
    out
end

@def vern9pre3 begin
    @vern9unpack
    b1Θdiff = @evalpoly(Θ, 6*r013, 24*r014, 60*r015, 120*r016, 210*r017, 336*r018, 504*r019)
    b8Θdiff = @evalpoly(Θ, 6*r083, 24*r084, 60*r085, 120*r086, 210*r087, 336*r088, 504*r089)
    b9Θdiff = @evalpoly(Θ, 6*r093, 24*r094, 60*r095, 120*r096, 210*r097, 336*r098, 504*r099)
    b10Θdiff = @evalpoly(Θ, 6*r103, 24*r104, 60*r105, 120*r106, 210*r107, 336*r108,
                         504*r109)
    b11Θdiff = @evalpoly(Θ, 6*r113, 24*r114, 60*r115, 120*r116, 210*r117, 336*r118,
                         504*r119)
    b12Θdiff = @evalpoly(Θ, 6*r123, 24*r124, 60*r125, 120*r126, 210*r127, 336*r128,
                         504*r129)
    b13Θdiff = @evalpoly(Θ, 6*r133, 24*r134, 60*r135, 120*r136, 210*r137, 336*r138,
                         504*r139)
    b14Θdiff = @evalpoly(Θ, 6*r143, 24*r144, 60*r145, 120*r146, 210*r147, 336*r148,
                         504*r149)
    b15Θdiff = @evalpoly(Θ, 6*r153, 24*r154, 60*r155, 120*r156, 210*r157, 336*r158,
                         504*r159)
    b17Θdiff = @evalpoly(Θ, 6*r173, 24*r174, 60*r175, 120*r176, 210*r177, 336*r178,
                         504*r179)
    b18Θdiff = @evalpoly(Θ, 6*r183, 24*r184, 60*r185, 120*r186, 210*r187, 336*r188,
                         504*r189)
    b19Θdiff = @evalpoly(Θ, 6*r193, 24*r194, 60*r195, 120*r196, 210*r197, 336*r198,
                         504*r199)
    b20Θdiff = @evalpoly(Θ, 6*r203, 24*r204, 60*r205, 120*r206, 210*r207, 336*r208,
                         504*r209)
    b21Θdiff = @evalpoly(Θ, 6*r213, 24*r214, 60*r215, 120*r216, 210*r217, 336*r218,
                         504*r219)
    b22Θdiff = @evalpoly(Θ, 6*r223, 24*r224, 60*r225, 120*r226, 210*r227, 336*r228,
                         504*r229)
    b23Θdiff = @evalpoly(Θ, 6*r233, 24*r234, 60*r235, 120*r236, 210*r237, 336*r238,
                         504*r239)
    b24Θdiff = @evalpoly(Θ, 6*r243, 24*r244, 60*r245, 120*r246, 210*r247, 336*r248,
                         504*r249)
    b25Θdiff = @evalpoly(Θ, 6*r253, 24*r254, 60*r255, 120*r256, 210*r257, 336*r258,
                         504*r259)
    b26Θdiff = @evalpoly(Θ, 6*r263, 24*r264, 60*r265, 120*r266, 210*r267, 336*r268,
                         504*r269)
    invdt2 = inv(dt)^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache},
                                  idxs::Nothing, T::Type{Val{3}})
    @vern9pre3
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds (k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
                      k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff +
                      k[8] * b14Θdiff +
                      k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
                      k[13] * b19Θdiff +
                      k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
                      k[17] * b23Θdiff +
                      k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff) * invdt2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                  T::Type{Val{3}})
    @vern9pre3
    return (k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
            k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
            k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
            k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
            k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
            k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
            k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
            k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
            k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff) * invdt2
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache},
                                   idxs::Nothing, T::Type{Val{3}})
    @vern9pre3
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
                                       k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
                                       k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
                                       k[11] * b17Θdiff + k[12] * b18Θdiff +
                                       k[13] * b19Θdiff + k[14] * b20Θdiff +
                                       k[15] * b21Θdiff + k[16] * b22Θdiff +
                                       k[17] * b23Θdiff + k[18] * b24Θdiff +
                                       k[19] * b25Θdiff + k[20] * b26Θdiff) * invdt2
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                   T::Type{Val{3}})
    @vern9pre3
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
                                    k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
                                    k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
                                    k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
                                    k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
                                    k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
                                    k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
                                    k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
                                    k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
                                    k[20][idxs] * b26Θdiff) * invdt2
    out
end

@def vern9pre4 begin
    @vern9unpack
    b1Θdiff = @evalpoly(Θ, 24*r014, 120*r015, 360*r016, 840*r017, 1680*r018, 3024*r019)
    b8Θdiff = @evalpoly(Θ, 24*r084, 120*r085, 360*r086, 840*r087, 1680*r088, 3024*r089)
    b9Θdiff = @evalpoly(Θ, 24*r094, 120*r095, 360*r096, 840*r097, 1680*r098, 3024*r099)
    b10Θdiff = @evalpoly(Θ, 24*r104, 120*r105, 360*r106, 840*r107, 1680*r108, 3024*r109)
    b11Θdiff = @evalpoly(Θ, 24*r114, 120*r115, 360*r116, 840*r117, 1680*r118, 3024*r119)
    b12Θdiff = @evalpoly(Θ, 24*r124, 120*r125, 360*r126, 840*r127, 1680*r128, 3024*r129)
    b13Θdiff = @evalpoly(Θ, 24*r134, 120*r135, 360*r136, 840*r137, 1680*r138, 3024*r139)
    b14Θdiff = @evalpoly(Θ, 24*r144, 120*r145, 360*r146, 840*r147, 1680*r148, 3024*r149)
    b15Θdiff = @evalpoly(Θ, 24*r154, 120*r155, 360*r156, 840*r157, 1680*r158, 3024*r159)
    b17Θdiff = @evalpoly(Θ, 24*r174, 120*r175, 360*r176, 840*r177, 1680*r178, 3024*r179)
    b18Θdiff = @evalpoly(Θ, 24*r184, 120*r185, 360*r186, 840*r187, 1680*r188, 3024*r189)
    b19Θdiff = @evalpoly(Θ, 24*r194, 120*r195, 360*r196, 840*r197, 1680*r198, 3024*r199)
    b20Θdiff = @evalpoly(Θ, 24*r204, 120*r205, 360*r206, 840*r207, 1680*r208, 3024*r209)
    b21Θdiff = @evalpoly(Θ, 24*r214, 120*r215, 360*r216, 840*r217, 1680*r218, 3024*r219)
    b22Θdiff = @evalpoly(Θ, 24*r224, 120*r225, 360*r226, 840*r227, 1680*r228, 3024*r229)
    b23Θdiff = @evalpoly(Θ, 24*r234, 120*r235, 360*r236, 840*r237, 1680*r238, 3024*r239)
    b24Θdiff = @evalpoly(Θ, 24*r244, 120*r245, 360*r246, 840*r247, 1680*r248, 3024*r249)
    b25Θdiff = @evalpoly(Θ, 24*r254, 120*r255, 360*r256, 840*r257, 1680*r258, 3024*r259)
    b26Θdiff = @evalpoly(Θ, 24*r264, 120*r265, 360*r266, 840*r267, 1680*r268, 3024*r269)
    invdt3 = inv(dt)^3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache},
                                  idxs::Nothing, T::Type{Val{4}})
    @vern9pre4
    #@.. broadcast=false k[1]*b1Θdiff + k[2]*b8Θdiff + k[3]*b9Θdiff + k[4]*b10Θdiff + k[5]*b11Θdiff + k[6]*b12Θdiff + k[7]*b13Θdiff + k[8]*b14Θdiff + k[9]*b15Θdiff + k[11]*b17Θdiff + k[12]*b18Θdiff + k[13]*b19Θdiff + k[14]*b20Θdiff + k[15]*b21Θdiff + k[16]*b22Θdiff + k[17]*b23Θdiff + k[18]*b24Θdiff + k[19]*b25Θdiff + k[20]*b26Θdiff
    return @inbounds (k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff + k[4] * b10Θdiff +
                      k[5] * b11Θdiff + k[6] * b12Θdiff + k[7] * b13Θdiff +
                      k[8] * b14Θdiff +
                      k[9] * b15Θdiff + k[11] * b17Θdiff + k[12] * b18Θdiff +
                      k[13] * b19Θdiff +
                      k[14] * b20Θdiff + k[15] * b21Θdiff + k[16] * b22Θdiff +
                      k[17] * b23Θdiff +
                      k[18] * b24Θdiff + k[19] * b25Θdiff + k[20] * b26Θdiff) * invdt3
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                  T::Type{Val{4}})
    @vern9pre4
    return (k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff + k[3][idxs] * b9Θdiff +
            k[4][idxs] * b10Θdiff + k[5][idxs] * b11Θdiff +
            k[6][idxs] * b12Θdiff + k[7][idxs] * b13Θdiff +
            k[8][idxs] * b14Θdiff + k[9][idxs] * b15Θdiff +
            k[11][idxs] * b17Θdiff + k[12][idxs] * b18Θdiff +
            k[13][idxs] * b19Θdiff + k[14][idxs] * b20Θdiff +
            k[15][idxs] * b21Θdiff + k[16][idxs] * b22Θdiff +
            k[17][idxs] * b23Θdiff + k[18][idxs] * b24Θdiff +
            k[19][idxs] * b25Θdiff + k[20][idxs] * b26Θdiff) * invdt3
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache},
                                   idxs::Nothing, T::Type{Val{4}})
    @vern9pre4
    @inbounds @.. broadcast=false out=(k[1] * b1Θdiff + k[2] * b8Θdiff + k[3] * b9Θdiff +
                                       k[4] * b10Θdiff + k[5] * b11Θdiff + k[6] * b12Θdiff +
                                       k[7] * b13Θdiff + k[8] * b14Θdiff + k[9] * b15Θdiff +
                                       k[11] * b17Θdiff + k[12] * b18Θdiff +
                                       k[13] * b19Θdiff + k[14] * b20Θdiff +
                                       k[15] * b21Θdiff + k[16] * b22Θdiff +
                                       k[17] * b23Θdiff + k[18] * b24Θdiff +
                                       k[19] * b25Θdiff + k[20] * b26Θdiff) * invdt3
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{Vern9ConstantCache, Vern9Cache}, idxs,
                                   T::Type{Val{4}})
    @vern9pre4
    @views @.. broadcast=false out=(k[1][idxs] * b1Θdiff + k[2][idxs] * b8Θdiff +
                                    k[3][idxs] * b9Θdiff + k[4][idxs] * b10Θdiff +
                                    k[5][idxs] * b11Θdiff + k[6][idxs] * b12Θdiff +
                                    k[7][idxs] * b13Θdiff + k[8][idxs] * b14Θdiff +
                                    k[9][idxs] * b15Θdiff + k[11][idxs] * b17Θdiff +
                                    k[12][idxs] * b18Θdiff + k[13][idxs] * b19Θdiff +
                                    k[14][idxs] * b20Θdiff + k[15][idxs] * b21Θdiff +
                                    k[16][idxs] * b22Θdiff + k[17][idxs] * b23Θdiff +
                                    k[18][idxs] * b24Θdiff + k[19][idxs] * b25Θdiff +
                                    k[20][idxs] * b26Θdiff) * invdt3
    out
end

"""

"""
@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
                                  T::Type{Val{0}})
    Θ1 = 1 - Θ
    # return @.. broadcast=false y₀ + dt*Θ*(k[1] + Θ1*(k[2] + Θ*(k[3]+Θ1*(k[4] + Θ*(k[5] + Θ1*(k[6]+Θ*k[7]))))))
    return @inbounds y₀ +
                     dt * Θ *
                     (k[1] +
                      Θ1 * (k[2] +
                       Θ * (k[3] + Θ1 * (k[4] + Θ * (k[5] + Θ1 * (k[6] + Θ * k[7]))))))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP8ConstantCache, DP8Cache}, idxs,
                                  T::Type{Val{0}})
    Θ1 = 1 - Θ
    # return @.. broadcast=false y₀[idxs] + dt*Θ*(k[1][idxs] + Θ1*(k[2][idxs] + Θ*(k[3][idxs]+Θ1*(k[4][idxs] + Θ*(k[5][idxs] + Θ1*(k[6][idxs]+Θ*k[7][idxs]))))))
    return y₀[idxs] +
           dt * Θ *
           (k[1][idxs] +
            Θ1 * (k[2][idxs] +
             Θ * (k[3][idxs] +
              Θ1 * (k[4][idxs] + Θ * (k[5][idxs] + Θ1 * (k[6][idxs] + Θ * k[7][idxs]))))))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
                                   T::Type{Val{0}})
    Θ1 = 1 - Θ
    @inbounds @.. broadcast=false out=y₀ +
                                      dt * Θ *
                                      (k[1] +
                                       Θ1 * (k[2] +
                                        Θ * (k[3] +
                                         Θ1 * (k[4] + Θ * (k[5] + Θ1 * (k[6] + Θ * k[7]))))))
    #@inbounds for i in eachindex(out)
    #  out[i] = y₀[i] + dt*Θ*(k[1][i] + Θ1*(k[2][i] + Θ*(k[3][i]+Θ1*(k[4][i] + Θ*(k[5][i] + Θ1*(k[6][i]+Θ*k[7][i]))))))
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP8ConstantCache, DP8Cache}, idxs,
                                   T::Type{Val{0}})
    Θ1 = 1 - Θ
    @views @.. broadcast=false out=y₀[idxs] +
                                   dt * Θ *
                                   (k[1][idxs] +
                                    Θ1 * (k[2][idxs] +
                                     Θ * (k[3][idxs] +
                                      Θ1 * (k[4][idxs] +
                                       Θ *
                                       (k[5][idxs] + Θ1 * (k[6][idxs] + Θ * k[7][idxs]))))))
    #@inbounds for (j,i) in enumerate(idxs)
    #  out[j] = y₀[i] + dt*Θ*(k[1][i] + Θ1*(k[2][i] + Θ*(k[3][i]+Θ1*(k[4][i] + Θ*(k[5][i] + Θ1*(k[6][i]+Θ*k[7][i]))))))
    #end
    out
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
                                  T::Type{Val{1}})
    @inbounds b1diff = @.. broadcast=false k[1]+k[2]
    @inbounds b2diff = @.. broadcast=false -2*k[2]+2*k[3]+2*k[4]
    @inbounds b3diff = @.. broadcast=false -3 * k[3]-6 * k[4]+3*k[5]+3*k[6]
    @inbounds b4diff = @.. broadcast=false 4 * k[4] - 8 * k[5] - 12 * k[6]+4 * k[7]
    @inbounds b5diff = @.. broadcast=false 5 * k[5] + 15 * k[6]-15 * k[7]
    @inbounds b6diff = @.. broadcast=false -6 * k[6]+18 * k[7]
    @inbounds b7diff = @.. broadcast=false -7*k[7]
    # return @.. broadcast=false b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff + Θ*(b5diff + Θ*(b6diff + Θ*b7diff)))))
    return b1diff +
           Θ *
           (b2diff + Θ * (b3diff + Θ * (b4diff + Θ * (b5diff + Θ * (b6diff + Θ * b7diff)))))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DP8ConstantCache, DP8Cache}, idxs,
                                  T::Type{Val{1}})
    b1diff = @.. broadcast=false k[1][idxs]+k[2][idxs]
    b2diff = @.. broadcast=false -2*k[2][idxs]+2*k[3][idxs]+2*k[4][idxs]
    b3diff = @.. broadcast=false -3 * k[3][idxs]-6 * k[4][idxs]+3*k[5][idxs]+3*k[6][idxs]
    b4diff = @.. broadcast=false 4 * k[4][idxs] - 8 * k[5][idxs] -
                                 12 * k[6][idxs]+4 * k[7][idxs]
    b5diff = @.. broadcast=false 5 * k[5][idxs] + 15 * k[6][idxs]-15 * k[7][idxs]
    b6diff = @.. broadcast=false -6 * k[6][idxs]+18 * k[7][idxs]
    b7diff = @.. broadcast=false -7*k[7][idxs]
    # return @.. broadcast=false b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff + Θ*(b5diff + Θ*(b6diff + Θ*b7diff)))))
    return b1diff +
           Θ *
           (b2diff + Θ * (b3diff + Θ * (b4diff + Θ * (b5diff + Θ * (b6diff + Θ * b7diff)))))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP8ConstantCache, DP8Cache}, idxs::Nothing,
                                   T::Type{Val{1}})
    # b1diff = k[1] + k[2]
    # b2diff = -2*k[2] + 2*k[3] + 2*k[4]
    # b3diff = -3*k[3] - 6*k[4] + 3*k[5] + 3*k[6]
    # b4diff = 4*k[4] - 8*k[5] - 12*k[6] + 4*k[7]
    # b5diff = 5*k[5] + 15*k[6] - 15*k[7]
    # b6diff = -6*k[6] + 18*k[7]
    # @.. broadcast=false out = b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff +
    #                                                Θ*(b5diff + Θ*(b6diff - 7*k[7]*Θ)))))
    @views @.. broadcast=false out=k[1] + k[2] +
                                   Θ * (-2 * k[2] + 2 * k[3] + 2 * k[4] +
                                    Θ * (-3 * k[3] - 6 * k[4] + 3 * k[5] + 3 * k[6] +
                                     Θ * (4 * k[4] - 8 * k[5] - 12 * k[6] + 4 * k[7] +
                                      Θ * (5 * k[5] + 15 * k[6] - 15 * k[7] +
                                       Θ * (-6 * k[6] + 18 * k[7] - 7 * k[7] * Θ)))))
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DP8ConstantCache, DP8Cache}, idxs,
                                   T::Type{Val{1}})
    # b1diff = k[1][idxs] + k[2][idxs]
    # b2diff = -2*k[2][idxs] + 2*k[3][idxs] + 2*k[4][idxs]
    # b3diff = -3*k[3][idxs] - 6*k[4][idxs] + 3*k[5][idxs] + 3*k[6][idxs]
    # b4diff = 4*k[4][idxs] - 8*k[5][idxs] - 12*k[6][idxs] + 4*k[7][idxs]
    # b5diff = 5*k[5][idxs] + 15*k[6][idxs] - 15*k[7][idxs]
    # b6diff = -6*k[6][idxs] + 18*k[7][idxs]
    #@views @.. broadcast=false out = b1diff + Θ*(b2diff + Θ*(b3diff + Θ*(b4diff +
    #                                               Θ*(b5diff + Θ*(b6diff - 7*k[7][idxs]*Θ)))))
    @views @.. broadcast=false out=k[1][idxs] + k[2][idxs] +
                                   Θ * (-2 * k[2][idxs] + 2 * k[3][idxs] + 2 * k[4][idxs] +
                                    Θ *
                                    (-3 * k[3][idxs] - 6 * k[4][idxs] + 3 * k[5][idxs] +
                                     3 * k[6][idxs] +
                                     Θ *
                                     (4 * k[4][idxs] - 8 * k[5][idxs] - 12 * k[6][idxs] +
                                      4 * k[7][idxs] +
                                      Θ *
                                      (5 * k[5][idxs] + 15 * k[6][idxs] - 15 * k[7][idxs] +
                                       Θ * (-6 * k[6][idxs] + 18 * k[7][idxs] -
                                        7 * k[7][idxs] * Θ)))))
    out
end

"""

"""
@def dprkn6unpack begin if typeof(cache) <: OrdinaryDiffEqMutableCache
    @unpack r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41, r54, r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10, rp34, rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51, rp64, rp63, rp62, rp61 = cache.tab
else
    @unpack r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41, r54, r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10, rp34, rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51, rp64, rp63, rp62, rp61 = cache
end end

@def dprkn6pre0 begin
    @dprkn6unpack
    b1Θ = @evalpoly(Θ, r10, r11, r12, r13, r14)
    b3Θ = @evalpoly(Θ, 0, r31, r32, r33, r34)
    b4Θ = @evalpoly(Θ, 0, r41, r42, r43, r44)
    b5Θ = @evalpoly(Θ, 0, r51, r52, r53, r54)
    b6Θ = @evalpoly(Θ, 0, r61, r62, r63, r64)

    bp1Θ = @evalpoly(Θ, rp10, rp11, rp12, rp13, rp14)
    bp3Θ = @evalpoly(Θ, 0, rp31, rp32, rp33, rp34)
    bp4Θ = @evalpoly(Θ, 0, rp41, rp42, rp43, rp44)
    bp5Θ = @evalpoly(Θ, 0, rp51, rp52, rp53, rp54)
    bp6Θ = @evalpoly(Θ, 0, rp61, rp62, rp63, rp64)

    kk1, kk2, kk3 = k
    k1, k2 = kk1.x
    k3, k4 = kk2.x
    k5, k6 = kk3.x

    duprev, uprev = y₀.x
    dtsq = dt^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DPRKN6ConstantCache, DPRKN6Cache},
                                  idxs::Nothing, T::Type{Val{0}})
    @dprkn6pre0
    return ArrayPartition(duprev +
                          dt * Θ *
                          (bp1Θ * k1 + bp3Θ * k3 +
                           bp4Θ * k4 + bp5Θ * k5 + bp6Θ * k6),
                          uprev +
                          dt * Θ *
                          (duprev +
                           dt * Θ * (b1Θ * k1 + b3Θ * k3 +
                                     b4Θ * k4 + b5Θ * k5 + b6Θ * k6)))
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
                                  cache::Union{DPRKN6ConstantCache, DPRKN6Cache}, idxs,
                                  T::Type{Val{0}})
    @dprkn6pre0
    return ArrayPartition(duprev[idxs] +
                          dt * Θ *
                          (bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
                           bp4Θ * k4[idxs] + bp5Θ * k5[idxs] + bp6Θ * k6[idxs]),
                          uprev[idxs] +
                          dt * Θ *
                          (duprev[idxs] +
                           dt * Θ *
                           (b1Θ * k1[idxs] +
                            b3Θ * k3[idxs] +
                            b4Θ * k4[idxs] + b5Θ * k5[idxs] + b6Θ * k6[idxs])))
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DPRKN6ConstantCache, DPRKN6Cache},
                                   idxs::Nothing, T::Type{Val{0}})
    @dprkn6pre0
    @inbounds @.. broadcast=false out.x[2]=uprev +
                                           dt * Θ *
                                           (duprev +
                                            dt * Θ *
                                            (b1Θ * k1 +
                                             b3Θ * k3 +
                                             b4Θ * k4 + b5Θ * k5 + b6Θ * k6))
    @inbounds @.. broadcast=false out.x[1]=duprev +
                                           dt * Θ *
                                           (bp1Θ * k1 + bp3Θ * k3 +
                                            bp4Θ * k4 + bp5Θ * k5 + bp6Θ * k6)
    #for i in eachindex(out.x[1])
    #  out.x[2][i]  = uprev[i] + dt*Θ*(duprev[i] + dt*Θ*(b1Θ*k1[i] +
    #                                                    b3Θ*k3[i] +
    #                                                    b4Θ*k4[i] + b5Θ*k5[i] + b6Θ*k6[i]))
    #  out.x[1][i] =  duprev[i] + dt*Θ*(bp1Θ*k1[i] + bp3Θ*k3[i] +
    #                                   bp4Θ*k4[i] + bp5Θ*k5[i] + bp6Θ*k6[i])
    #end
    out
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
                                   cache::Union{DPRKN6ConstantCache, DPRKN6Cache}, idxs,
                                   T::Type{Val{0}})
    @dprkn6pre0
    @inbounds @.. broadcast=false out.x[2]=uprev[idxs] +
                                           dt * Θ *
                                           (duprev[idxs] +
                                            dt * Θ *
                                            (b1Θ * k1[idxs] +
                                             b3Θ * k3[idxs] +
                                             b4Θ * k4[idxs] + b5Θ * k5[idxs] +
                                             b6Θ * k6[idxs]))
    @inbounds @.. broadcast=false out.x[1]=duprev[idxs] +
                                           dt * Θ *
                                           (bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
                                            bp4Θ * k4[idxs] + bp5Θ * k5[idxs] +
                                            bp6Θ * k6[idxs])
    #for (j,i) in enumerate(idxs)
    #  out.x[2][j]  = uprev[i] + dt*Θ*(duprev[i] + dt*Θ*(b1Θ*k1[i] +
    #                                                    b3Θ*k3[i] +
    #                                                    b4Θ*k4[i] + b5Θ*k5[i] + b6Θ*k6[i]))
    #  out.x[1][j] = duprev[i] + dt*Θ*(bp1Θ*k1[i] + bp3Θ*k3[i] +
    #                                  bp4Θ*k4[i] + bp5Θ*k5[i] + bp6Θ*k6[i])
    #end
    out
end
