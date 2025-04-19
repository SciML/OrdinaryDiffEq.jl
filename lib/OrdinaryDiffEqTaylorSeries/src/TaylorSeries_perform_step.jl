using TaylorDiff: TaylorDiff, extract_derivative, extract_derivative!

@inline make_taylor(all::Vararg{X, P}) where {P, X <: AbstractArray} = TaylorArray(
    Base.first(all), Base.tail(all))
@inline make_taylor(all::Vararg{X, P}) where {P, X} = TaylorScalar(all)

function initialize!(integrator, cache::ExplicitTaylor2ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylor2ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    k1 = f(uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    k2 = f(u1, p, t1).partials[1]
    u = @.. uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylor2Cache)
    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k1, k2, k3, utilde, tmp = cache

    # The following code is written to be fully non-allocating
    f(k1, uprev, p, t)
    u1 = make_taylor(uprev, k1)
    t1 = TaylorScalar{1}(t, one(t))
    out1 = make_taylor(k1, k2)
    f(out1, u1, p, t1)
    @.. u = uprev + dt * k1 + dt^2 / 2 * k2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end

function initialize!(integrator, cache::ExplicitTaylorConstantCache{P}) where {P}
    integrator.kshortsize = P
    integrator.k = typeof(integrator.k)(undef, P)
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylorConstantCache{P}, repeat_step = false) where {P}
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack jet = cache
    utaylor = jet(uprev, t)
    u = map(x -> evaluate_polynomial(x, dt), utaylor)
    if integrator.opts.adaptive
        utilde = TaylorDiff.get_coefficient(utaylor, P) * dt^P
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P + 1)
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylorCache{P}) where {P}
    integrator.kshortsize = P
    resize!(integrator.k, P)
    # Setup k pointers
    for i in 1:P
        integrator.k[i] = get_coefficient(cache.utaylor, i)
    end
    return nothing
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylorCache{P}, repeat_step = false) where {P}
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack jet, utaylor, utilde, tmp, atmp, thread = cache

    jet(utaylor, uprev, t)
    for i in eachindex(utaylor)
        u[i] = @inline evaluate_polynomial(utaylor[i], dt)
    end
    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=TaylorDiff.get_coefficient(utaylor, P) *
                                                 dt^P
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P + 1)
    return nothing
end

function initialize!(integrator, cache::ExplicitTaylorAdaptiveOrderCache)
    integrator.kshortsize = cache.max_order
    resize!(integrator.k, cache.max_order)
    # Setup k pointers
    for i in 1:(cache.max_order)
        integrator.k[i] = get_coefficient(cache.utaylor, i)
    end
    return nothing
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylorAdaptiveOrderCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack jets, current_order, min_order, max_order, utaylor, utilde, tmp, atmp, thread = cache

    jet_index = current_order[] - min_order + 1
    # compute one additional order for adaptive order
    jet = jets[jet_index + 1]
    jet(utaylor, uprev, t)
    for i in eachindex(utaylor)
        u[i] = @inline evaluate_polynomial(utaylor[i], dt)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, current_order[] + 1)
    if integrator.opts.adaptive
        min_work = Inf
        start_order = max(min_order, current_order[] - 1)
        end_order = min(max_order, current_order[] + 1)
        for i in start_order:end_order
            A = i * i
            @.. broadcast=false thread=thread utilde=TaylorDiff.get_coefficient(
                utaylor, i) * dt^i
            calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            EEst = integrator.opts.internalnorm(atmp, t)

            # backup
            e = integrator.EEst
            qold = integrator.qold
            # calculate dt
            integrator.EEst = EEst
            dtpropose = step_accept_controller!(integrator, alg,
                stepsize_controller!(integrator, alg))
            # restore
            integrator.EEst = e
            integrator.qold = qold

            work = A / dtpropose
            if work < min_work
                cache.current_order[] = i
                min_work = work
                integrator.EEst = EEst
            end
        end
    end
    return nothing
end

function initialize!(integrator, cache::ExplicitTaylorAdaptiveOrderConstantCache)
    integrator.kshortsize = cache.max_order
    integrator.k = typeof(integrator.k)(undef, cache.max_order)
    return nothing
end

@muladd function perform_step!(
        integrator, cache::ExplicitTaylorAdaptiveOrderConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, false)
    @unpack jets, current_order, min_order, max_order = cache

    jet_index = current_order[] - min_order + 1
    # compute one additional order for adaptive order
    jet = jets[jet_index + 1]
    utaylor = jet(uprev, t)
    u = map(x -> evaluate_polynomial(x, dt), utaylor)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, current_order[] + 1)
    if integrator.opts.adaptive
        min_work = Inf
        start_order = max(min_order, current_order[] - 1)
        end_order = min(max_order, current_order[] + 1)
        for i in start_order:end_order
            A = i * i
            utilde = TaylorDiff.get_coefficient(utaylor, i) * dt^i
            atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            EEst = integrator.opts.internalnorm(atmp, t)

            # backup
            e = integrator.EEst
            qold = integrator.qold
            # calculate dt
            integrator.EEst = EEst
            dtpropose = step_accept_controller!(integrator, alg,
                stepsize_controller!(integrator, alg))
            # restore
            integrator.EEst = e
            integrator.qold = qold

            work = A / dtpropose
            if work < min_work
                cache.current_order[] = i
                min_work = work
                integrator.EEst = EEst
            end
        end
    end
    return nothing
end

function initialize!(integrator, cache::ImplicitTaylor2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ImplicitTaylor2ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
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

    if integrator.opts.adaptive && integrator.differential_vars !== nothing
        atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                  integrator.opts.abstol
        integrator.EEst += integrator.opts.internalnorm(atmp, t)
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end
