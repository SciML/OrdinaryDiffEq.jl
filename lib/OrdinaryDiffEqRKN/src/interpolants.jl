@def dprkn6unpack begin
    if cache isa OrdinaryDiffEqMutableCache
        @unpack r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41, r54, r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10, rp34, rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51, rp64, rp63, rp62, rp61 = cache.tab
    else
        @unpack r14, r13, r12, r11, r10, r34, r33, r32, r31, r44, r43, r42, r41, r54, r53, r52, r51, r64, r63, r62, r61, rp14, rp13, rp12, rp11, rp10, rp34, rp33, rp32, rp31, rp44, rp43, rp42, rp41, rp54, rp53, rp52, rp51, rp64, rp63, rp62, rp61 = cache
    end
end

@def dprkn6pre0 begin
    @dprkn6unpack
    b1Θ = @evalpoly(Θ, r10, r11, r12, r13, r14)
    b3Θ = Θ * @evalpoly(Θ, r31, r32, r33, r34)
    b4Θ = Θ * @evalpoly(Θ, r41, r42, r43, r44)
    b5Θ = Θ * @evalpoly(Θ, r51, r52, r53, r54)
    b6Θ = Θ * @evalpoly(Θ, r61, r62, r63, r64)

    bp1Θ = @evalpoly(Θ, rp10, rp11, rp12, rp13, rp14)
    bp3Θ = Θ * @evalpoly(Θ, rp31, rp32, rp33, rp34)
    bp4Θ = Θ * @evalpoly(Θ, rp41, rp42, rp43, rp44)
    bp5Θ = Θ * @evalpoly(Θ, rp51, rp52, rp53, rp54)
    bp6Θ = Θ * @evalpoly(Θ, rp61, rp62, rp63, rp64)

    kk1, kk2, kk3 = k
    k1, k2 = kk1.x
    k3, k4 = kk2.x
    k5, k6 = kk3.x

    duprev, uprev = y₀.x
    dtsq = dt^2
end

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DPRKN6ConstantCache, DPRKN6Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
    @dprkn6pre0
    return ArrayPartition(
        duprev +
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
        T::Type{Val{0}}, differential_vars::Nothing)
    @dprkn6pre0
    return ArrayPartition(
        duprev[idxs] +
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

@muladd function _ode_interpolant(Θ, dt, y₀, y₁, k,
        cache::Union{DPRKN6ConstantCache, DPRKN6Cache}, idxs::Number,
        T::Type{Val{0}}, differential_vars::Nothing)
    @dprkn6pre0
    halfsize = length(y₀) ÷ 2
    if idxs <= halfsize
        duprev[idxs] +
        dt * Θ *
        (bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
         bp4Θ * k4[idxs] + bp5Θ * k5[idxs] + bp6Θ * k6[idxs])
    else
        idxs = idxs - halfsize
        uprev[idxs] +
        dt * Θ *
        (duprev[idxs] +
         dt * Θ *
         (b1Θ * k1[idxs] +
          b3Θ * k3[idxs] +
          b4Θ * k4[idxs] + b5Θ * k5[idxs] + b6Θ * k6[idxs]))
    end
end

@muladd function _ode_interpolant!(out, Θ, dt, y₀, y₁, k,
        cache::Union{DPRKN6ConstantCache, DPRKN6Cache},
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing)
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
        T::Type{Val{0}}, differential_vars::Nothing)
    @dprkn6pre0
    halfsize = length(y₀) ÷ 2
    isfirsthalf = idxs .<= halfsize
    secondhalf = idxs .> halfsize
    firstidxs = idxs[isfirsthalf]
    secondidxs_shifted = idxs[secondhalf]
    secondidxs = secondidxs_shifted .- halfsize

    @views @.. broadcast=false out[secondhalf]=uprev[secondidxs] +
                                               dt * Θ *
                                               (duprev[secondidxs] +
                                                dt * Θ *
                                                (b1Θ * k1[secondidxs] +
                                                 b3Θ * k3[secondidxs] +
                                                 b4Θ * k4[secondidxs] +
                                                 b5Θ * k5[secondidxs] +
                                                 b6Θ * k6[secondidxs]))
    @views @.. broadcast=false out[isfirsthalf]=duprev[firstidxs] +
                                                dt * Θ *
                                                (bp1Θ * k1[firstidxs] +
                                                 bp3Θ * k3[firstidxs] +
                                                 bp4Θ * k4[firstidxs] +
                                                 bp5Θ * k5[firstidxs] +
                                                 bp6Θ * k6[firstidxs])
    out
end
