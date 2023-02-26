function initialize!(integrator,
                     cache::Union{ImplicitEulerConstantCache,
                                  ImplicitMidpointConstantCache,
                                  TrapezoidConstantCache,
                                  TRBDF2ConstantCache,
                                  SDIRK2ConstantCache,
                                  SDIRK22ConstantCache,
                                  SSPSDIRK2ConstantCache,
                                  Cash4ConstantCache,
                                  Hairer4ConstantCache,
                                  ESDIRK54I8L2SAConstantCache,
                                  ESDIRK436L2SA2ConstantCache,
                                  ESDIRK437L2SAConstantCache,
                                  ESDIRK547L2SA2ConstantCache,
                                  SFSDIRK4ConstantCache,
                                  SFSDIRK5ConstantCache,
                                  SFSDIRK6ConstantCache,
                                  SFSDIRK7ConstantCache,
                                  SFSDIRK8ConstantCache})
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.destats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator,
                     cache::Union{ImplicitEulerCache,
                                  ImplicitMidpointCache,
                                  TrapezoidCache,
                                  TRBDF2Cache,
                                  SDIRK2Cache,
                                  SDIRK22Cache,
                                  SSPSDIRK2Cache,
                                  Cash4Cache,
                                  Hairer4Cache,
                                  ESDIRK54I8L2SACache,
                                  ESDIRK436L2SA2Cache,
                                  ESDIRK437L2SACache,
                                  ESDIRK547L2SA2Cache,
                                  SFSDIRK4Cache,
                                  SFSDIRK5Cache,
                                  SFSDIRK6Cache,
                                  SFSDIRK7Cache,
                                  SFSDIRK8Cache})
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ImplicitEulerConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else # :constant
        nlsolver.z = zero(u)
    end

    nlsolver.tmp = uprev
    nlsolver.γ = 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + z

    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        tmp = r *
              integrator.opts.internalnorm.((u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
        atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    else
        integrator.EEst = 1
    end

    integrator.fsallast = f(u, p, t + dt)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ImplicitEulerCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;atmp, nlsolver) = cache
    (;z, tmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast=false z=dt * integrator.fsalfirst
    else # :constant
        z .= zero(eltype(u))
    end

    nlsolver.tmp .= uprev
    nlsolver.γ = 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=uprev + z

    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        @.. broadcast=false tmp=r * integrator.opts.internalnorm((u - uprev) / dt1 -
                                                             (uprev - uprev2) / dt2, t)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    else
        integrator.EEst = 1
    end
    integrator.destats.nf += 1
    f(integrator.fsallast, u, p, t + dt)
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    γ = 1 // 2
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else # :constant
        nlsolver.z = zero(u)
    end

    nlsolver.tmp = uprev
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + z

    integrator.fsallast = f(u, p, t + dt)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ImplicitMidpointCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;nlsolver) = cache
    (;z, tmp) = nlsolver
    mass_matrix = integrator.f.mass_matrix
    alg = unwrap_alg(integrator, true)
    γ = 1 // 2
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast=false z=dt * integrator.fsalfirst
    else # :constant
        z .= zero(eltype(u))
    end

    nlsolver.tmp = uprev
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=nlsolver.tmp + z

    integrator.destats.nf += 1
    f(integrator.fsallast, u, p, t + dt)
end

@muladd function perform_step!(integrator, cache::TrapezoidConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    # precalculations
    γ = 1 // 2
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess: constant extrapolation
    nlsolver.z = uprev

    if f.mass_matrix === I
        nlsolver.tmp = @.. broadcast=false uprev * inv(γdt)+integrator.fsalfirst
    else
        nlsolver.tmp = (f.mass_matrix * uprev) .* inv(γdt) .+ integrator.fsalfirst
    end
    nlsolver.α = 1
    nlsolver.γ = γ
    nlsolver.method = COEFFICIENT_MULTISTEP
    u = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    if integrator.opts.adaptive
        if integrator.iter > 2
            # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
            # use 3rd divided differences (DD) a la SPICE and Shampine

            # TODO: check numerical stability
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2

            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
            r = c * dt^3 / 2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

            # tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
            DD31 = (u - uprev) / dt1 - (uprev - uprev2) / dt2
            DD30 = (uprev - uprev2) / dt3 - (uprev2 - uprev3) / dt4
            tmp = r * integrator.opts.internalnorm((DD31 - DD30) / dt5, t)
            atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                                       integrator.opts.reltol, integrator.opts.internalnorm,
                                       t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
            if integrator.EEst <= 1
                cache.uprev3 = uprev2
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            integrator.EEst = 1
            cache.uprev3 = integrator.uprev2
            cache.tprev2 = integrator.tprev
        else
            integrator.EEst = 1
        end
    end

    integrator.fsallast = f(u, p, t + dt)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::TrapezoidCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;atmp, nlsolver) = cache
    (;z, tmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    γ = 1 // 2
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess: constant extrapolation
    @.. broadcast=false z=uprev
    invγdt = inv(γdt)
    if mass_matrix === I
        @.. broadcast=false tmp=uprev * invγdt + integrator.fsalfirst
    else
        mul!(u, mass_matrix, uprev)
        @.. broadcast=false tmp=u * invγdt + integrator.fsalfirst
    end
    nlsolver.α = 1
    nlsolver.γ = γ
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=z

    if integrator.opts.adaptive
        if integrator.iter > 2
            # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
            # use 3rd divided differences (DD) a la SPICE and Shampine

            # TODO: check numerical stability
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2

            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
            r = c * dt^3 / 2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

            # @.. broadcast=false tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
            @.. broadcast=false tmp=r * integrator.opts.internalnorm((((u - uprev) / dt1 -
                                                                   (uprev - uprev2) / dt2) #DD31
                                                                  -
                                                                  ((uprev - uprev2) / dt3 -
                                                                   (uprev2 - uprev3) / dt4)) /
                                                                 dt5, t)
            calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
            if integrator.EEst <= 1
                copyto!(cache.uprev3, uprev2)
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            integrator.EEst = 1
            copyto!(cache.uprev3, integrator.uprev2)
            cache.tprev2 = integrator.tprev
        else
            integrator.EEst = 1
        end
    end

    integrator.destats.nf += 1
    f(integrator.fsallast, u, p, t + dt)
end

@muladd function perform_step!(integrator, cache::TRBDF2ConstantCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, d, ω, btilde1, btilde2, btilde3, α1, α2) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # FSAL
    zprev = dt * integrator.fsalfirst

    ##### Solve Trapezoid Step

    # TODO: Add extrapolation
    zᵧ = zprev
    nlsolver.z = zᵧ
    nlsolver.c = γ

    nlsolver.tmp = uprev + d * zprev
    zᵧ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve BDF2 Step

    ### Initial Guess From Shampine
    z = α1 * zprev + α2 * zᵧ
    nlsolver.z = z
    nlsolver.c = 1

    nlsolver.tmp = uprev + ω * zprev + ω * zᵧ
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + d * z

    ################################### Finalize

    if integrator.opts.adaptive
        tmp = btilde1 * zprev + btilde2 * zᵧ + btilde3 * z
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.destats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::TRBDF2Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;zprev, zᵧ, atmp, nlsolver) = cache
    (;z, tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    b = nlsolver.ztmp
    (;γ, d, ω, btilde1, btilde2, btilde3, α1, α2) = cache.tab
    alg = unwrap_alg(integrator, true)

    # FSAL
    @.. broadcast=false zprev=dt * integrator.fsalfirst
    markfirststage!(nlsolver)

    ##### Solve Trapezoid Step

    # TODO: Add extrapolation
    @.. broadcast=false zᵧ=zprev
    z .= zᵧ
    @.. broadcast=false tmp=uprev + d * zprev
    nlsolver.c = γ
    zᵧ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve BDF2 Step

    ### Initial Guess From Shampine
    @.. broadcast=false z=α1 * zprev + α2 * zᵧ
    @.. broadcast=false tmp=uprev + ω * zprev + ω * zᵧ
    nlsolver.c = 1
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + d * z

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * zprev + btilde2 * zᵧ + btilde3 * z
        if alg.smooth_est && isnewton(nlsolver) # From Shampine
            est = nlsolver.cache.dz
            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                                linu = _vec(est))

            integrator.destats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z / dt
end

@muladd function perform_step!(integrator, cache::TRBDF2Cache{<:Array}, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;zprev, zᵧ, atmp, nlsolver) = cache
    (;z, tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    b = nlsolver.ztmp
    (;γ, d, ω, btilde1, btilde2, btilde3, α1, α2) = cache.tab
    alg = unwrap_alg(integrator, true)

    # FSAL
    @inbounds @simd ivdep for i in eachindex(u)
        zprev[i] = dt * integrator.fsalfirst[i]
    end
    markfirststage!(nlsolver)

    ##### Solve Trapezoid Step

    # TODO: Add extrapolation
    copyto!(zᵧ, zprev)
    copyto!(z, zᵧ)
    @inbounds @simd ivdep for i in eachindex(u)
        tmp[i] = uprev[i] + d * zprev[i]
    end
    nlsolver.c = γ
    zᵧ .= nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve BDF2 Step

    ### Initial Guess From Shampine
    @inbounds @simd ivdep for i in eachindex(u)
        z[i] = α1 * zprev[i] + α2 * zᵧ[i]
    end
    @inbounds @simd ivdep for i in eachindex(u)
        tmp[i] = uprev[i] + ω * zprev[i] + ω * zᵧ[i]
    end
    nlsolver.c = 1
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = tmp[i] + d * z[i]
    end

    ################################### Finalize

    if integrator.opts.adaptive
        @inbounds @simd ivdep for i in eachindex(u)
            tmp[i] = btilde1 * zprev[i] + btilde2 * zᵧ[i] + btilde3 * z[i]
        end
        if alg.smooth_est && isnewton(nlsolver) # From Shampine
            est = nlsolver.cache.dz

            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                                linu = _vec(est))

            integrator.destats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @inbounds @simd ivdep for i in eachindex(u)
        integrator.fsallast[i] = z[i] / dt
    end
end

@muladd function perform_step!(integrator, cache::SDIRK2ConstantCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if integrator.success_iter > 0 && !integrator.reeval_fsal &&
       alg.extrapolant == :interpolant
        current_extrapolant!(u, t + dt, integrator)
        z₁ = u - uprev
    elseif alg.extrapolant == :linear
        z₁ = dt * integrator.fsalfirst
    else
        z₁ = zero(u)
    end

    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ### Initial Guess Is α₁ = c₂/γ, c₂ = 0 => z₂ = α₁z₁ = 0
    z₂ = zero(u)
    nlsolver.z = z₂
    nlsolver.tmp = uprev - z₁
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = uprev + z₁ / 2 + z₂ / 2

    ################################### Finalize

    if integrator.opts.adaptive
        tmp = z₁ / 2 - z₂ / 2
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.destats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = f(u, p, t)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SDIRK2Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if integrator.success_iter > 0 && !integrator.reeval_fsal &&
       alg.extrapolant == :interpolant
        current_extrapolant!(u, t + dt, integrator)
        @.. broadcast=false z₁=u - uprev
    elseif alg.extrapolant == :linear
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    else
        z₁ .= zero(eltype(u))
    end
    nlsolver.z = z₁

    ##### Step 1
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 2

    ### Initial Guess Is α₁ = c₂/γ, c₂ = 0 => z₂ = α₁z₁ = 0
    z₂ .= zero(eltype(u))
    nlsolver.z = z₂
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    @.. broadcast=false tmp=uprev - z₁
    nlsolver.tmp = tmp
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=uprev + z₁ / 2 + z₂ / 2

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=z₁ / 2 - z₂ / 2
        if alg.smooth_est && isnewton(nlsolver) # From Shampine
            est = nlsolver.cache.dz
            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                                linu = _vec(est))
            integrator.destats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.destats.nf += 1
    f(integrator.fsallast, u, p, t)
end

@muladd function perform_step!(integrator, cache::SDIRK22ConstantCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;a, α, β) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    # precalculations
    γ = a * dt
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess
    zprev = dt * integrator.fsalfirst
    nlsolver.z = zprev

    # first stage
    nlsolver.tmp = uprev + γdt * integrator.fsalfirst
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    uprev = α * nlsolver.tmp + β * z

    # final stage
    γ = dt
    γdt = γ * dt
    markfirststage!(nlsolver)
    nlsolver.tmp = uprev + γdt * integrator.fsalfirst
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp

    if integrator.opts.adaptive
        if integrator.iter > 2
            # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
            # use 3rd divided differences (DD) a la SPICE and Shampine

            # TODO: check numerical stability
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2

            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
            r = c * dt^3 / 2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

            DD31 = (u - uprev) / dt1 - (uprev - uprev2) / dt2
            DD30 = (uprev - uprev2) / dt3 - (uprev2 - uprev3) / dt4
            tmp = r * abs((DD31 - DD30) / dt5)
            atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
                                       integrator.opts.reltol, integrator.opts.internalnorm,
                                       t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
            if integrator.EEst <= 1
                cache.uprev3 = uprev2
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            integrator.EEst = 1
            cache.uprev3 = integrator.uprev2
            cache.tprev2 = integrator.tprev
        else
            integrator.EEst = 1
        end
    end

    integrator.destats.nf += 2
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SDIRK22Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;atmp, nlsolver) = cache
    (;z, tmp) = nlsolver
    (;a, α, β) = cache.tab
    alg = unwrap_alg(integrator, true)
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    γ = a * dt
    γdt = γ * dt
    markfirststage!(nlsolver)

    # first stage
    @.. broadcast=false z=dt * integrator.fsalfirst
    @.. broadcast=false tmp=uprev + γdt * integrator.fsalfirst
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=α * tmp + β * z

    # final stage
    γ = dt
    γdt = γ * dt
    markfirststage!(nlsolver)
    @.. broadcast=false tmp=uprev + γdt * integrator.fsalfirst
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=nlsolver.tmp

    if integrator.opts.adaptive
        if integrator.iter > 2
            # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
            # use 3rd divided differences (DD) a la SPICE and Shampine

            # TODO: check numerical stability
            uprev2 = integrator.uprev2
            tprev = integrator.tprev
            uprev3 = cache.uprev3
            tprev2 = cache.tprev2

            dt1 = dt * (t + dt - tprev)
            dt2 = (t - tprev) * (t + dt - tprev)
            dt3 = (t - tprev) * (t - tprev2)
            dt4 = (tprev - tprev2) * (t - tprev2)
            dt5 = t + dt - tprev2
            c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
            r = c * dt^3 / 2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

            @inbounds for i in eachindex(u)
                DD31 = (u[i] - uprev[i]) / dt1 - (uprev[i] - uprev2[i]) / dt2
                DD30 = (uprev[i] - uprev2[i]) / dt3 - (uprev2[i] - uprev3[i]) / dt4
                tmp[i] = r * abs((DD31 - DD30) / dt5)
            end
            calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
            if integrator.EEst <= 1
                copyto!(cache.uprev3, uprev2)
                cache.tprev2 = tprev
            end
        elseif integrator.success_iter > 0
            integrator.EEst = 1
            copyto!(cache.uprev3, integrator.uprev2)
            cache.tprev2 = integrator.tprev
        else
            integrator.EEst = 1
        end
    end

    integrator.destats.nf += 2
    f(integrator.fsallast, u, p, t + dt)
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    γ = eltype(u)(1 // 4)
    c2 = typeof(t)(3 // 4)

    markfirststage!(nlsolver)

    # initial guess
    if integrator.success_iter > 0 && !integrator.reeval_fsal &&
       alg.extrapolant == :interpolant
        current_extrapolant!(u, t + dt, integrator)
        z₁ = u - uprev
    elseif alg.extrapolant == :linear
        z₁ = dt * integrator.fsalfirst
    else
        z₁ = zero(u)
    end
    nlsolver.z = z₁

    ##### Step 1

    tstep = t + dt
    u = uprev + γ * z₁

    nlsolver.c = 1
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 2

    ### Initial Guess Is α₁ = c₂/γ
    z₂ = c2 / γ
    nlsolver.z = z₂

    nlsolver.tmp = uprev + z₁ / 2
    nlsolver.c = 1
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + z₂ / 2

    ################################### Finalize

    integrator.fsallast = f(u, p, t)
    integrator.destats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SSPSDIRK2Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, nlsolver) = cache
    (;tmp) = nlsolver
    alg = unwrap_alg(integrator, true)

    γ = eltype(u)(1 // 4)
    c2 = typeof(t)(3 // 4)
    markfirststage!(nlsolver)

    # initial guess
    if integrator.success_iter > 0 && !integrator.reeval_fsal &&
       alg.extrapolant == :interpolant
        current_extrapolant!(u, t + dt, integrator)
        @.. broadcast=false z₁=u - uprev
    elseif alg.extrapolant == :linear
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    else
        z₁ .= zero(eltype(u))
    end
    nlsolver.z = z₁
    nlsolver.tmp = uprev

    ##### Step 1
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 2

    ### Initial Guess Is α₁ = c₂/γ
    @.. broadcast=false z₂=c2 / γ
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + z₁ / 2
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + z₂ / 2

    ################################### Finalize

    integrator.destats.nf += 1
    f(integrator.fsallast, u, p, t)
end

@muladd function perform_step!(integrator, cache::Cash4ConstantCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c2, c3, c4) = cache.tab
    (;b1hat1, b2hat1, b3hat1, b4hat1, b1hat2, b2hat2, b3hat2, b4hat2) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z = z₁

    nlsolver.c = γ
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ = zero(u)
    nlsolver.z = z₂

    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    z₃ = z₁
    nlsolver.z = z₃

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    z₄ = z₃
    nlsolver.z = z₄

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use yhat2 for prediction
    z₅ = b1hat2 * z₁ + b2hat2 * z₂ + b3hat2 * z₃ + b4hat2 * z₄
    nlsolver.z = z₅

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = 1
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₅

    ################################### Finalize

    if integrator.opts.adaptive
        if alg.embedding == 3
            btilde1 = b1hat2 - a51
            btilde2 = b2hat2 - a52
            btilde3 = b3hat2 - a53
            btilde4 = b4hat2 - a54
            btilde5 = -γ
        else
            btilde1 = b1hat1 - a51
            btilde2 = b2hat1 - a52
            btilde3 = b3hat1 - a53
            btilde4 = b4hat1 - a54
            btilde5 = -γ
        end

        tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.destats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₅ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Cash4Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c2, c3, c4) = cache.tab
    (;b1hat1, b2hat1, b3hat1, b4hat1, b1hat2, b2hat2, b3hat2, b4hat2) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ .= zero(eltype(z₁))
    nlsolver.z = z₁
    nlsolver.c = γ
    nlsolver.tmp = uprev

    # initial step of NLNewton iteration
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(z₂))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    @.. broadcast=false z₃=z₁
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    @.. broadcast=false z₅=b1hat2 * z₁ + b2hat2 * z₂ + b3hat2 * z₃ + b4hat2 * z₄
    nlsolver.z = z₅
    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = 1
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₅

    ################################### Finalize

    if integrator.opts.adaptive
        if alg.embedding == 3
            btilde1 = b1hat2 - a51
            btilde2 = b2hat2 - a52
            btilde3 = b3hat2 - a53
            btilde4 = b4hat2 - a54
            btilde5 = -γ
        else
            btilde1 = b1hat1 - a51
            btilde2 = b2hat1 - a52
            btilde3 = b3hat1 - a53
            btilde4 = b4hat1 - a54
            btilde5 = -γ
        end

        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅
        if alg.smooth_est && isnewton(nlsolver) # From Shampine
            est = nlsolver.cache.dz
            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                                linu = _vec(est))
            integrator.destats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₅ / dt
end

@muladd function perform_step!(integrator, cache::SFSDIRK4ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c2, c3, c4) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z = z₁

    nlsolver.c = γ
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ = zero(u)
    nlsolver.z = z₂

    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    z₃ = z₁
    nlsolver.z = z₃

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    z₄ = z₃
    nlsolver.z = z₄

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Final Step

    u = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄

    integrator.fsallast = z₄ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SFSDIRK4Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c2, c3, c4) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)
    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ .= zero(eltype(z₁))
    nlsolver.z = z₁
    nlsolver.c = γ
    nlsolver.tmp = uprev

    # initial step of NLNewton iteration
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(z₂))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    @.. broadcast=false z₃=z₁
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄

    ################################### Finalize
    @.. broadcast=false integrator.fsallast=z₄ / dt
end

@muladd function perform_step!(integrator, cache::SFSDIRK5ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, c2, c3, c4, c5) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z = z₁

    nlsolver.c = γ
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ = zero(u)
    nlsolver.z = z₂

    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    z₃ = z₁
    nlsolver.z = z₃

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    z₄ = z₃
    nlsolver.z = z₄

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    z₅ = z₄
    nlsolver.z = z₅

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Final Step

    u = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅

    integrator.fsallast = z₅ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SFSDIRK5Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, c2, c3, c4, c5) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)
    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ .= zero(eltype(z₁))
    nlsolver.z = z₁
    nlsolver.c = γ
    nlsolver.tmp = uprev

    # initial step of NLNewton iteration
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(z₂))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    @.. broadcast=false z₃=z₁
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    @.. broadcast=false z₅=z₄
    nlsolver.z = z₅

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    ################################### Finalize
    @.. broadcast=false u=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅

    @.. broadcast=false integrator.fsallast=z₅ / dt
end

@muladd function perform_step!(integrator, cache::SFSDIRK6ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, c2, c3, c4, c5, c6) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z = z₁

    nlsolver.c = γ
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ = zero(u)
    nlsolver.z = z₂

    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    z₃ = z₁
    nlsolver.z = z₃

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    z₄ = z₃
    nlsolver.z = z₄

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    z₅ = z₄
    nlsolver.z = z₅

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    z₆ = z₅
    nlsolver.z = z₆

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Final Step

    u = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆

    integrator.fsallast = z₆ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SFSDIRK6Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, c2, c3, c4, c5, c6) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)
    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ .= zero(eltype(z₁))
    nlsolver.z = z₁
    nlsolver.c = γ
    nlsolver.tmp = uprev

    # initial step of NLNewton iteration
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(z₂))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    @.. broadcast=false z₃=z₁
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    @.. broadcast=false z₅=z₄
    nlsolver.z = z₅

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    # Use constant z prediction
    @.. broadcast=false z₆=z₅
    nlsolver.z = z₆

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################### Finalize
    @.. broadcast=false u=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                          a76 * z₆

    @.. broadcast=false integrator.fsallast=z₆ / dt
end

@muladd function perform_step!(integrator, cache::SFSDIRK7ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, c2, c3, c4, c5, c6, c7) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z = z₁

    nlsolver.c = γ
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ = zero(u)
    nlsolver.z = z₂

    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    z₃ = z₁
    nlsolver.z = z₃

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    z₄ = z₃
    nlsolver.z = z₄

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    z₅ = z₄
    nlsolver.z = z₅

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    # Use constant z prediction
    z₆ = z₅
    nlsolver.z = z₆

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    # Use constant z prediction
    z₇ = z₆
    nlsolver.z = z₇

    nlsolver.tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Final Step

    u = uprev + a81 * z₁ + a82 * z₂ + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ + a87 * z₇

    integrator.fsallast = z₇ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SFSDIRK7Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, z₇, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, c2, c3, c4, c5, c6, c7) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)
    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ .= zero(eltype(z₁))
    nlsolver.z = z₁
    nlsolver.c = γ
    nlsolver.tmp = uprev

    # initial step of NLNewton iteration
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(z₂))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    @.. broadcast=false z₃=z₁
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    @.. broadcast=false z₅=z₄
    nlsolver.z = z₅

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    # Use constant z prediction
    @.. broadcast=false z₆=z₅
    nlsolver.z = z₆

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    # Use constant z prediction
    @.. broadcast=false z₇=z₆
    nlsolver.z = z₇

    @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                            a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################### Finalize
    @.. broadcast=false u=uprev + a81 * z₁ + a82 * z₂ + a83 * z₃ + a84 * z₄ + a85 * z₅ +
                          a86 * z₆ + a87 * z₇

    @.. broadcast=false integrator.fsallast=z₇ / dt
end

@muladd function perform_step!(integrator, cache::SFSDIRK8ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, a91, a92, a93, a94, a95, a96, a97, a98, c2, c3, c4, c5, c6, c7, c8) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z = z₁

    nlsolver.c = γ
    nlsolver.tmp = uprev
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ = zero(u)
    nlsolver.z = z₂

    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    z₃ = z₁
    nlsolver.z = z₃

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    z₄ = z₃
    nlsolver.z = z₄

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    z₅ = z₄
    nlsolver.z = z₅

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    # Use constant z prediction
    z₆ = z₅
    nlsolver.z = z₆

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    # Use constant z prediction
    z₇ = z₆
    nlsolver.z = z₇

    nlsolver.tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    # Use constant z prediction
    z₈ = z₇
    nlsolver.z = z₈

    nlsolver.tmp = uprev + a81 * z₁ + a82 * z₂ + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ +
                   a87 * z₇
    nlsolver.c = c8
    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Final Step

    u = uprev + a91 * z₁ + a92 * z₂ + a93 * z₃ + a94 * z₄ + a95 * z₅ + a96 * z₆ + a97 * z₇ +
        a98 * z₈

    integrator.fsallast = z₈ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SFSDIRK8Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, z₇, z₈, nlsolver) = cache
    (;tmp) = nlsolver
    W = isnewton(nlsolver) ? get_W(nlsolver) : nothing
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, a91, a92, a93, a94, a95, a96, a97, a98, c2, c3, c4, c5, c6, c7, c8) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)
    ##### Step 1

    # TODO: Add extrapolation for guess
    z₁ .= zero(eltype(z₁))
    nlsolver.z = z₁
    nlsolver.c = γ
    nlsolver.tmp = uprev

    # initial step of NLNewton iteration
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(z₂))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    # Guess starts from z₁
    @.. broadcast=false z₃=z₁
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=z₃
    nlsolver.z = z₄

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use constant z prediction
    @.. broadcast=false z₅=z₄
    nlsolver.z = z₅

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    # Use constant z prediction
    @.. broadcast=false z₆=z₅
    nlsolver.z = z₆

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    # Use constant z prediction
    @.. broadcast=false z₇=z₆
    nlsolver.z = z₇

    @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                            a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    # Use constant z prediction
    @.. broadcast=false z₈=z₇
    nlsolver.z = z₈

    @.. broadcast=false tmp=uprev + a81 * z₁ + a82 * z₂ + a83 * z₃ + a84 * z₄ + a85 * z₅ +
                            a86 * z₆ + a87 * z₇
    nlsolver.c = c8
    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################### Finalize
    @.. broadcast=false u=uprev + a91 * z₁ + a92 * z₂ + a93 * z₃ + a94 * z₄ + a95 * z₅ +
                          a96 * z₆ + a97 * z₇ + a98 * z₈

    @.. broadcast=false integrator.fsallast=z₈ / dt
end

@muladd function perform_step!(integrator, cache::Hairer4ConstantCache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c2, c3, c4) = cache.tab
    (;α21, α31, α32, α41, α43) = cache.tab
    (;bhat1, bhat2, bhat3, bhat4, btilde1, btilde2, btilde3, btilde4, btilde5) = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # TODO: Add extrapolation for guess
    z₁ = zero(u)
    nlsolver.z, nlsolver.tmp = z₁, uprev
    nlsolver.c = γ
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    z₂ = α21 * z₁
    nlsolver.z = z₂
    nlsolver.tmp = uprev + a21 * z₁
    nlsolver.c = c2
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    z₃ = α31 * z₁ + α32 * z₂
    nlsolver.z = z₃
    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    z₄ = α41 * z₁ + α43 * z₃
    nlsolver.z = z₄
    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use yhat2 for prediction
    z₅ = bhat1 * z₁ + bhat2 * z₂ + bhat3 * z₃ + bhat4 * z₄
    nlsolver.z = z₅
    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = 1
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₅

    ################################### Finalize

    if integrator.opts.adaptive
        tmp = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅
        if isnewton(nlsolver) && alg.smooth_est # From Shampine
            integrator.destats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₅ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::Hairer4Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    (;γ, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, c2, c3, c4) = cache.tab
    (;α21, α31, α32, α41, α43) = cache.tab
    (;bhat1, bhat2, bhat3, bhat4, btilde1, btilde2, btilde3, btilde4, btilde5) = cache.tab
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if integrator.success_iter > 0 && !integrator.reeval_fsal &&
       alg.extrapolant == :interpolant
        current_extrapolant!(u, t + dt, integrator)
        @.. broadcast=false z₁=u - uprev
    elseif alg.extrapolant == :linear
        @.. broadcast=false z₁=dt * integrator.fsalfirst
    else
        z₁ .= zero(eltype(z₁))
    end
    nlsolver.z = z₁
    nlsolver.tmp = uprev

    ##### Step 1

    nlsolver.c = γ
    z₁ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ##### Step 2

    @.. broadcast=false z₂=α21 * z₁
    nlsolver.z = z₂
    @.. broadcast=false tmp=uprev + a21 * z₁
    nlsolver.tmp = tmp
    nlsolver.c = c2
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    @.. broadcast=false z₃=α31 * z₁ + α32 * z₂
    nlsolver.z = z₃
    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    @.. broadcast=false z₄=α41 * z₁ + α43 * z₃
    nlsolver.z = z₄
    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    # Use yhat prediction
    @.. broadcast=false z₅=bhat1 * z₁ + bhat2 * z₂ + bhat3 * z₃ + bhat4 * z₄
    nlsolver.z = z₅
    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = 1
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₅

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅
        if alg.smooth_est && isnewton(nlsolver) # From Shampine
            est = nlsolver.cache.dz
            linres = dolinsolve(integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                                linu = _vec(est))

            integrator.destats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(atmp, est, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₅ / dt
end

@muladd function perform_step!(integrator, cache::ESDIRK54I8L2SAConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    a71, a72, a73, a74, a75, a76,
    a81, a82, a83, a84, a85, a86, a87,
    c3, c4, c5, c6, c7,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, btilde8 = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # TODO: Add extrapolation for guess

    ##### Step 1

    z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = zero(z₁)

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.z = z₃ = zero(z₂)

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = zero(z₃)

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = z₅ = zero(z₄)

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = z₆ = zero(z₅)

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    nlsolver.z = z₇ = zero(z₆)

    nlsolver.tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    nlsolver.z = z₈ = zero(z₇)

    nlsolver.tmp = uprev + a81 * z₁ + a82 * z₂ + a83 * z₃ + a84 * z₄ + a85 * z₅ + a86 * z₆ +
                   a87 * z₇
    nlsolver.c = 1
    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₈

    ################################### Finalize

    if integrator.opts.adaptive
        est = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
              btilde6 * z₆ + btilde7 * z₇ + btilde8 * z₈
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₈ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK54I8L2SACache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, z₇, z₈, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    a71, a72, a73, a74, a75, a76,
    a81, a82, a83, a84, a85, a86, a87,
    c3, c4, c5, c6, c7,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, btilde8 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    ##### Step 1

    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(u))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    nlsolver.z = fill!(z₃, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    nlsolver.z = fill!(z₄, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = fill!(z₅, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = fill!(z₆, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    nlsolver.z = fill!(z₇, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                            a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 8

    nlsolver.z = fill!(z₈, zero(eltype(u)))

    @.. broadcast=false nlsolver.tmp=uprev + a81 * z₁ + a82 * z₂ + a83 * z₃ + a84 * z₄ +
                                     a85 * z₅ + a86 * z₆ + a87 * z₇
    nlsolver.c = oneunit(nlsolver.c)
    z₈ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₈

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇ + btilde8 * z₈
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₈ / dt
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK436L2SA2ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    c3, c4, c5, c6,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6 = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # TODO: Add extrapolation for guess

    ##### Step 1

    z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = zero(z₁)

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.z = z₃ = zero(z₂)

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = zero(z₃)

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = z₅ = zero(z₄)

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = z₆ = zero(z₅)

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₆

    ################################### Finalize

    if integrator.opts.adaptive
        est = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
              btilde6 * z₆
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₆ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK436L2SA2Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    c3, c4, c5, c6,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    ##### Step 1

    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(u))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    nlsolver.z = fill!(z₃, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    nlsolver.z = fill!(z₄, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = fill!(z₅, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = fill!(z₆, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₆

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅ + btilde6 * z₆
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₆ / dt
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK437L2SAConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    a71, a72, a73, a74, a75, a76,
    c3, c4, c5, c6, c7,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # TODO: Add extrapolation for guess

    ##### Step 1

    z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = zero(z₁)

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.z = z₃ = zero(z₂)

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = zero(z₃)

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = z₅ = zero(z₄)

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = z₆ = zero(z₅)

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    nlsolver.z = z₇ = zero(z₆)

    nlsolver.tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₇

    ################################### Finalize

    if integrator.opts.adaptive
        est = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
              btilde6 * z₆ + btilde7 * z₇
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₇ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK437L2SACache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, z₇, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    a71, a72, a73, a74, a75, a76,
    c3, c4, c5, c6, c7,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    ##### Step 1

    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(u))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    nlsolver.z = fill!(z₃, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    nlsolver.z = fill!(z₄, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = fill!(z₅, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = fill!(z₆, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = fill!(z₇, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                            a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₇

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₇ / dt
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK547L2SA2ConstantCache,
                               repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    a71, a72, a73, a74, a75, a76,
    c3, c4, c5, c6, c7,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # TODO: Add extrapolation for guess

    ##### Step 1

    z₁ = dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation choice
    nlsolver.z = z₂ = zero(z₁)

    nlsolver.tmp = uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 3

    nlsolver.z = z₃ = zero(z₂)

    nlsolver.tmp = uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    nlsolver.z = z₄ = zero(z₃)

    nlsolver.tmp = uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = z₅ = zero(z₄)

    nlsolver.tmp = uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = z₆ = zero(z₅)

    nlsolver.tmp = uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    nlsolver.z = z₇ = zero(z₆)

    nlsolver.tmp = uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ + a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = nlsolver.tmp + γ * z₇

    ################################### Finalize

    if integrator.opts.adaptive
        est = btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ + btilde5 * z₅ +
              btilde6 * z₆ + btilde7 * z₇
        atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = z₇ ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::ESDIRK547L2SA2Cache, repeat_step = false)
    (;t, dt, uprev, u, f, p) = integrator
    (;z₁, z₂, z₃, z₄, z₅, z₆, z₇, atmp, nlsolver) = cache
    (;tmp) = nlsolver
    @unpack γ,
    a31, a32,
    a41, a42, a43,
    a51, a52, a53, a54,
    a61, a62, a63, a64, a65,
    a71, a72, a73, a74, a75, a76,
    c3, c4, c5, c6, c7,
    btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7 = cache.tab
    alg = unwrap_alg(integrator, true)

    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    ##### Step 1

    @.. broadcast=false z₁=dt * integrator.fsalfirst

    ##### Step 2

    # TODO: Add extrapolation for guess
    z₂ .= zero(eltype(u))
    nlsolver.z = z₂

    @.. broadcast=false tmp=uprev + γ * z₁
    nlsolver.c = 2γ
    z₂ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    isnewton(nlsolver) && set_new_W!(nlsolver, false)

    ################################## Solve Step 3

    nlsolver.z = fill!(z₃, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a31 * z₁ + a32 * z₂
    nlsolver.c = c3
    z₃ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 4

    # Use constant z prediction
    nlsolver.z = fill!(z₄, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a41 * z₁ + a42 * z₂ + a43 * z₃
    nlsolver.c = c4
    z₄ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 5

    nlsolver.z = fill!(z₅, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a51 * z₁ + a52 * z₂ + a53 * z₃ + a54 * z₄
    nlsolver.c = c5
    z₅ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 6

    nlsolver.z = fill!(z₆, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a61 * z₁ + a62 * z₂ + a63 * z₃ + a64 * z₄ + a65 * z₅
    nlsolver.c = c6
    z₆ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    ################################## Solve Step 7

    nlsolver.z = fill!(z₇, zero(eltype(u)))

    @.. broadcast=false tmp=uprev + a71 * z₁ + a72 * z₂ + a73 * z₃ + a74 * z₄ + a75 * z₅ +
                            a76 * z₆
    nlsolver.c = c7
    z₇ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false u=tmp + γ * z₇

    ################################### Finalize

    if integrator.opts.adaptive
        @.. broadcast=false tmp=btilde1 * z₁ + btilde2 * z₂ + btilde3 * z₃ + btilde4 * z₄ +
                                btilde5 * z₅ + btilde6 * z₆ + btilde7 * z₇
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
                             integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=z₇ / dt
    return
end
