@muladd function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::Tsit5Cache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false)

    if length(k) < 7 || always_calc_begin
        T = constvalue(recursive_unitless_bottom_eltype(u))
        T2 = constvalue(typeof(one(t)))
        @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
        @unpack k, tmp = cache
        
        # Loop-based implementation
        for i in 2:7
            @.. tmp = uprev + dt * sum(a[i, j] * k[j] for j in 1:(i-1))
            f(k[i], tmp, p, t + c[i-1] * dt)
        end
        
        # Copy results back into `k`
        for i in 1:7
            copyat_or_push!(k, i, k[i])
        end
    end
    nothing
end


@muladd function _ode_addsteps!(k, t, uprev, u, dt, f, p, cache::Tsit5ConstantCache,
                                always_calc_begin = false, allow_calc_end = true,
                                force_calc_end = false)
    if length(k) < 7 || always_calc_begin
        T = constvalue(recursive_unitless_bottom_eltype(u))
        T2 = constvalue(typeof(one(t)))
        @OnDemandTableauExtract Tsit5ConstantCacheActual T T2

        a_coeffs = [a[2, 1],
                    a[3, 1], a[3, 2],
                    a[4, 1], a[4, 2], a[4, 3],
                    a[5, 1], a[5, 2], a[5, 3], a[5, 4],
                    a[6, 1], a[6, 2], a[6, 3], a[6, 4], a[6, 5],
                    a[7, 1], a[7, 2], a[7, 3], a[7, 4], a[7, 5], a[7, 6]]

        c_coeffs = [c[1], c[2], c[3], c[4]]

        for i in 1:6
            u_increment = uprev
            for j in 1:(i-1)
                u_increment += dt * a_coeffs[(i*(i-1))÷2 + j] * k[j]
            end

            u_tmp = uprev + dt * u_increment
            t_tmp = t + (i <= 4 ? c_coeffs[i] : dt)
            copyat_or_push!(k, i, f(u_tmp, p, t_tmp))
        end

        utmp = uprev +
               dt *
               (a[7, 1] * k[1] + a[7, 2] * k[2] + a[7, 3] * k[3] + a[7, 4] * k[4] + a[7, 5] * k[5] + a[7, 6] * k[6])
        copyat_or_push!(k, 7, f(utmp, p, t + dt))
    end
    nothing
end

#=
@muladd function _ode_addsteps!(k,t,uprev,u,dt,f,p,cache::Tsit5Cache,always_calc_begin = false,allow_calc_end = true,force_calc_end = false)
  if length(k)<7 || always_calc_begin
    @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76 = cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,tmp = cache
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
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Tsit5ConstantCacheActual T T2

    k = Vector{typeof(k1)}(undef, 7)
    k[1] = integrator.fsalfirst

    for i in 2:6
        u_increment = uprev
        for j in 1:(i-1)
            u_increment += dt * a[i, j] * k[j]
        end
        k[i] = f(u_increment, p, t + c[i-1] * dt)
    end

    g6 = uprev
    for j in 1:5
        g6 += dt * a[6, j] * k[j]
    end
    k[6] = f(g6, p, t + dt)

    u = uprev
    for j in 1:6
        u += dt * a[7, j] * k[j]
    end

    integrator.fsallast = f(u, p, t + dt)
    k[7] = integrator.fsallast

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)

    if integrator.alg isa CompositeAlgorithm
        g7 = u
        integrator.eigen_est = integrator.opts.internalnorm(
            maximum(abs.((k[7] .- k[6]) ./ (g7 .- g6))), t)
    end

    if integrator.opts.adaptive
        utilde = dt * sum(btilde[j] * k[j] for j in 1:7)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    for i in 1:7
        integrator.k[i] = k[i]
    end
    integrator.u = u
end

function initialize!(integrator, cache::Tsit5Cache)
    integrator.kshortsize = 7
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.k[1]
    integrator.k[2] = cache.k[2]
    integrator.k[3] = cache.k[3]
    integrator.k[4] = cache.k[4]
    integrator.k[5] = cache.k[5]
    integrator.k[6] = cache.k[6]
    integrator.k[7] = cache.k[7]
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function perform_step!(integrator, cache::Tsit5Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
    @unpack k, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread = cache

    # Loop through stages 2 to 6
    for i in 2:6
        tmp .= uprev
        for j in 1:(i-1)
            tmp .+= dt * a[i, j] * k[j]
        end
        stage_limiter!(tmp, f, p, t + c[i-1] * dt)
        f(k[i], tmp, p, t + c[i-1] * dt)
    end

    # Handle the final (7th) stage separately
    u .= uprev
    for j in 1:6
        u .+= dt * a[7, j] * k[j]
    end
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k[7], u, p, t + dt)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)

    if integrator.alg isa CompositeAlgorithm
        g7 = u
        g6 = tmp
        @.. broadcast=false thread=thread utilde .= abs.((k[7] .- k[6]) ./ (g7 .- g6))
        integrator.eigen_est = integrator.opts.internalnorm(norm(utilde, Inf), t)
    end

    if integrator.opts.adaptive
        utilde .= 0.0
        for j in 1:7
            utilde .+= dt * btilde[j] * k[j]
        end
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t,
                             thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    return nothing
end