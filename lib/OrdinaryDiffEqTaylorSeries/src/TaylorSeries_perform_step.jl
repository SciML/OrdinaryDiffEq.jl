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
    (; t, dt, uprev, u, f, p) = integrator
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
    (; t, dt, uprev, u, f, p) = integrator
    (; k1, k2, k3, utilde, tmp) = cache

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
    (; t, dt, uprev, u, f, p) = integrator
    (; jet) = cache
    utaylor = jet(uprev, t)
    u = map(x -> evaluate_polynomial(x, dt), utaylor)
    if integrator.opts.adaptive
        utilde = TaylorDiff.get_coefficient(utaylor, P) * dt^(P + 1)
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
    (; t, dt, uprev, u, f, p) = integrator
    (; jet, utaylor, utilde, tmp, atmp, thread) = cache

    jet(utaylor, uprev, t)
    for i in eachindex(utaylor)
        u[i] = @inline evaluate_polynomial(utaylor[i], dt)
    end
    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=TaylorDiff.get_coefficient(utaylor, P) *
                                                 dt^(P + 1)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P + 1)
    return nothing
end
