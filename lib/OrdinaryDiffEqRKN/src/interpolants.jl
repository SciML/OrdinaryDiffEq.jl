@def dprkn6pre0 begin
    R, Rp = cache.tab.R, cache.tab.Rp
    b1Θ = @evalpoly(Θ, R[1, 1], R[1, 2], R[1, 3], R[1, 4], R[1, 5])
    b3Θ = @evalpoly(Θ, R[3, 1], R[3, 2], R[3, 3], R[3, 4], R[3, 5])
    b4Θ = @evalpoly(Θ, R[4, 1], R[4, 2], R[4, 3], R[4, 4], R[4, 5])
    b5Θ = @evalpoly(Θ, R[5, 1], R[5, 2], R[5, 3], R[5, 4], R[5, 5])
    b6Θ = @evalpoly(Θ, R[6, 1], R[6, 2], R[6, 3], R[6, 4], R[6, 5])

    bp1Θ = @evalpoly(Θ, Rp[1, 1], Rp[1, 2], Rp[1, 3], Rp[1, 4], Rp[1, 5])
    bp3Θ = @evalpoly(Θ, Rp[3, 1], Rp[3, 2], Rp[3, 3], Rp[3, 4], Rp[3, 5])
    bp4Θ = @evalpoly(Θ, Rp[4, 1], Rp[4, 2], Rp[4, 3], Rp[4, 4], Rp[4, 5])
    bp5Θ = @evalpoly(Θ, Rp[5, 1], Rp[5, 2], Rp[5, 3], Rp[5, 4], Rp[5, 5])
    bp6Θ = @evalpoly(Θ, Rp[6, 1], Rp[6, 2], Rp[6, 3], Rp[6, 4], Rp[6, 5])

    kk1, kk2, kk3 = k
    k1, k2 = kk1.x
    k3, k4 = kk2.x
    k5, k6 = kk3.x

    duprev, uprev = y₀.x
    dtsq = dt^2
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::DPRKN6Caches,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @dprkn6pre0
    return ArrayPartition(
        duprev +
            dt * Θ *
            (
            bp1Θ * k1 + bp3Θ * k3 +
                bp4Θ * k4 + bp5Θ * k5 + bp6Θ * k6
        ),
        uprev +
            dt * Θ *
            (
            duprev +
                dt * Θ * (
                b1Θ * k1 + b3Θ * k3 +
                    b4Θ * k4 + b5Θ * k5 + b6Θ * k6
            )
        )
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::DPRKN6Caches, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @dprkn6pre0
    return ArrayPartition(
        duprev[idxs] +
            dt * Θ *
            (
            bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
                bp4Θ * k4[idxs] + bp5Θ * k5[idxs] + bp6Θ * k6[idxs]
        ),
        uprev[idxs] +
            dt * Θ *
            (
            duprev[idxs] +
                dt * Θ *
                (
                b1Θ * k1[idxs] +
                    b3Θ * k3[idxs] +
                    b4Θ * k4[idxs] + b5Θ * k5[idxs] + b6Θ * k6[idxs]
            )
        )
    )
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::DPRKN6Caches, idxs::Number,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @dprkn6pre0
    halfsize = length(y₀) ÷ 2
    if idxs <= halfsize
        duprev[idxs] +
            dt * Θ *
            (
            bp1Θ * k1[idxs] + bp3Θ * k3[idxs] +
                bp4Θ * k4[idxs] + bp5Θ * k5[idxs] + bp6Θ * k6[idxs]
        )
    else
        idxs = idxs - halfsize
        uprev[idxs] +
            dt * Θ *
            (
            duprev[idxs] +
                dt * Θ *
                (
                b1Θ * k1[idxs] +
                    b3Θ * k3[idxs] +
                    b4Θ * k4[idxs] + b5Θ * k5[idxs] + b6Θ * k6[idxs]
            )
        )
    end
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::DPRKN6Caches,
        idxs::Nothing, T::Type{Val{0}}, differential_vars::Nothing
    )
    @dprkn6pre0
    @inbounds @.. broadcast = false out.x[2] = uprev +
        dt * Θ *
        (
        duprev +
            dt * Θ *
            (
            b1Θ * k1 +
                b3Θ * k3 +
                b4Θ * k4 + b5Θ * k5 + b6Θ * k6
        )
    )
    @inbounds @.. broadcast = false out.x[1] = duprev +
        dt * Θ *
        (
        bp1Θ * k1 + bp3Θ * k3 +
            bp4Θ * k4 + bp5Θ * k5 + bp6Θ * k6
    )
    #for i in eachindex(out.x[1])
    #  out.x[2][i]  = uprev[i] + dt*Θ*(duprev[i] + dt*Θ*(b1Θ*k1[i] +
    #                                                    b3Θ*k3[i] +
    #                                                    b4Θ*k4[i] + b5Θ*k5[i] + b6Θ*k6[i]))
    #  out.x[1][i] =  duprev[i] + dt*Θ*(bp1Θ*k1[i] + bp3Θ*k3[i] +
    #                                   bp4Θ*k4[i] + bp5Θ*k5[i] + bp6Θ*k6[i])
    #end
    out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::DPRKN6Caches, idxs,
        T::Type{Val{0}}, differential_vars::Nothing
    )
    @dprkn6pre0
    halfsize = length(y₀) ÷ 2
    isfirsthalf = idxs .<= halfsize
    secondhalf = idxs .> halfsize
    firstidxs = idxs[isfirsthalf]
    secondidxs_shifted = idxs[secondhalf]
    secondidxs = secondidxs_shifted .- halfsize

    @views @.. broadcast = false out[secondhalf] = uprev[secondidxs] +
        dt * Θ *
        (
        duprev[secondidxs] +
            dt * Θ *
            (
            b1Θ * k1[secondidxs] +
                b3Θ * k3[secondidxs] +
                b4Θ * k4[secondidxs] +
                b5Θ * k5[secondidxs] +
                b6Θ * k6[secondidxs]
        )
    )
    @views @.. broadcast = false out[isfirsthalf] = duprev[firstidxs] +
        dt * Θ *
        (
        bp1Θ * k1[firstidxs] +
            bp3Θ * k3[firstidxs] +
            bp4Θ * k4[firstidxs] +
            bp5Θ * k5[firstidxs] +
            bp6Θ * k6[firstidxs]
    )
    out
end
