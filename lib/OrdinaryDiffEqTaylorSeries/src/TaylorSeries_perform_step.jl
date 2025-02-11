using TaylorDiff: TaylorDiff, extract_derivative, extract_derivative!

@inline make_taylor(all::Vararg{X, P}) where {P, X <: AbstractArray} = TaylorArray(Base.first(all), Base.tail(all))
@inline make_taylor(all::Vararg{X, P}) where {P, X} = TaylorScalar(all)

function initialize!(integrator, cache::ExplicitTaylor2ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2ConstantCache, repeat_step = false)
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

function initialize!(integrator, cache::ExplicitTaylorConstantCache{P}) where P
    integrator.kshortsize = P
    integrator.k = typeof(integrator.k)(undef, P)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylorConstantCache{P}, repeat_step = false) where P
    @unpack t, dt, uprev, u, f, p = integrator
    us = typeof(u)[]
    integrator.k[1] = f(uprev, p, t)
    push!(us, integrator.k[1])
    u = @.. uprev + dt * us[1]
    dti = dt
    for i in 1:P-1
        ui = make_taylor(uprev, us...)
        ti = TaylorScalar{i}(t, one(t))
        integrator.k[i + 1] = f(ui, p, ti).partials[i]
        push!(us, integrator.k[i + 1] / (i + 1))
        dti *= dt
        u += dti * us[i + 1]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, P)
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitTaylorCache{P}) where P
    integrator.kshortsize = P
    resize!(integrator.k, P)
    # Setup k pointers
    for (i, k) in enumerate(cache.ks)
        integrator.k[i] = k
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::ExplicitTaylorCache{P}, repeat_step = false) where P
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack ks, us = cache

    # The following code is written to be fully non-allocating
    f(ks[1], uprev, p, t)
    @.. us[1] .= ks[1]
    @.. u = uprev + dt * us[1]
    dti = dt
    for i in 1:P-1
        ui = make_taylor(uprev, us[1:i]...)
        ti = TaylorScalar{i}(t, one(t))
        outi = make_taylor(ks[1:i+1]...)
        f(outi, ui, p, ti)
        us[i + 1] .= ks[i + 1] / (i + 1)
        dti *= dt
        @.. u += dti * us[i + 1]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end
