@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Tsit5CacheType,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 7 || always_calc_begin
        T = constvalue(recursive_unitless_bottom_eltype(u))
        T2 = constvalue(typeof(one(t)))
        @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
        (; k1, k2, k3, k4, k5, k6, k7, tmp) = cache
        @.. broadcast = false tmp = uprev + dt * (a21 * k1)
        f(k2, tmp, p, t + c1 * dt)
        @.. broadcast = false tmp = uprev + dt * (a31 * k1 + a32 * k2)
        f(k3, tmp, p, t + c2 * dt)
        @.. broadcast = false tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
        f(k4, tmp, p, t + c3 * dt)
        @.. broadcast = false tmp = uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
        f(k5, tmp, p, t + c4 * dt)
        @.. broadcast = false tmp = uprev +
            dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
        f(k6, tmp, p, t + dt)
        @.. broadcast = false tmp = uprev +
            dt *
            (
            a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 +
                a76 * k6
        )
        f(k7, tmp, p, t + dt)
        copyat_or_push!(k, 1, k1)
        copyat_or_push!(k, 2, k2)
        copyat_or_push!(k, 3, k3)
        copyat_or_push!(k, 4, k4)
        copyat_or_push!(k, 5, k5)
        copyat_or_push!(k, 6, k6)
        copyat_or_push!(k, 7, k7)
    end
    nothing
end

@muladd function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::Tsit5ConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 7 || always_calc_begin
        T = constvalue(recursive_unitless_bottom_eltype(u))
        T2 = constvalue(typeof(one(t)))
        @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
        copyat_or_push!(k, 1, f(uprev, p, t))
        copyat_or_push!(k, 2, f(uprev + dt * (a21 * k[1]), p, t + c1 * dt))
        copyat_or_push!(k, 3, f(uprev + dt * (a31 * k[1] + a32 * k[2]), p, t + c2 * dt))
        copyat_or_push!(
            k, 4,
            f(
                uprev + dt * (a41 * k[1] + a42 * k[2] + a43 * k[3]), p,
                t + c3 * dt
            )
        )
        copyat_or_push!(
            k, 5,
            f(
                uprev + dt * (a51 * k[1] + a52 * k[2] + a53 * k[3] + a54 * k[4]),
                p, t + c4 * dt
            )
        )
        copyat_or_push!(
            k, 6,
            f(
                uprev +
                    dt *
                    (a61 * k[1] + a62 * k[2] + a63 * k[3] + a64 * k[4] + a65 * k[5]),
                p, t + dt
            )
        )
        utmp = uprev +
            dt *
            (a71 * k[1] + a72 * k[2] + a73 * k[3] + a74 * k[4] + a75 * k[5] + a76 * k[6])
        copyat_or_push!(k, 7, f(utmp, p, t + dt))
    end
    nothing
end

#=
@muladd function _ode_addsteps!(k,t,uprev,u,dt,f,p,cache::Tsit5CacheType,always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<7 || always_calc_begin
    (; c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76) = cache.tab
    (; k1,k2,k3,k4,k5,k6,k7,tmp) = cache
    uidx = eachindex(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a21*k1[i])
    end
    f(k2,tmp,p,t+c1*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    f(k3,tmp,p,t+c2*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    f(k4,tmp,p,t+c3*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    f(k5,tmp,p,t+c4*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    f(k6,tmp,p,t+dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    f(k7,u,p,t+dt)
    copyat_or_push!(k,1,k1)
    copyat_or_push!(k,2,k2)
    copyat_or_push!(k,3,k3)
    copyat_or_push!(k,4,k4)
    copyat_or_push!(k,5,k5)
    copyat_or_push!(k,6,k6)
    copyat_or_push!(k,7,k7)
  end
  nothing
end
=#

function initialize!(integrator, cache::Tsit5ConstantCache)
    integrator.kshortsize = 7
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    g6 = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(g6, p, t + dt)
    u = uprev + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    integrator.fsallast = f(u, p, t + dt)
    k7 = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    if integrator.alg isa CompositeAlgorithm
        g7 = u
        # Hairer II, page 22 modified to use the Inf norm
        integrator.eigen_est = integrator.opts.internalnorm(
            maximum(abs.((k7 .- k6) ./ (g7 .- g6))), t
        )
    end
    if integrator.opts.adaptive
        utilde = dt *
            (
            btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 +
                btilde6 * k6 + btilde7 * k7
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.u = u
end

function initialize!(integrator, cache::Tsit5CacheType)
    integrator.kshortsize = 7
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function perform_step!(integrator, cache::Tsit5CacheType, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
    (; k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    a = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a * k1
    stage_limiter!(tmp, f, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, f, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, f, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, f, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, f, p, t + dt)
    f(k6, tmp, p, t + dt)
    @.. broadcast = false thread = thread u = uprev +
        dt * (
        a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 +
            a75 * k5 + a76 * k6
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k7, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    if integrator.alg isa CompositeAlgorithm
        g7 = u
        g6 = tmp
        # Hairer II, page 22 modified to use Inf norm
        @.. broadcast = false thread = thread utilde = abs((k7 - k6) / (g7 - g6))
        integrator.eigen_est = integrator.opts.internalnorm(norm(utilde, Inf), t)
    end
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde2 * k2 +
                btilde3 * k3 + btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    return nothing
end
