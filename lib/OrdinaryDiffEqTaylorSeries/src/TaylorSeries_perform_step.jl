
function initialize!(integrator, cache::ExplicitTaylor2ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end

@muladd function perform_step!(integrator, cache::ExplicitTaylor2ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    k1 = f(uprev, p, t)
    k2 = TaylorDiff.derivative(_u -> f(_u, p, t), uprev, k1, Val(1))
    k3 = TaylorDiff.derivative(_t -> f(uprev, p, _t), t, Val(1))
    u = @.. uprev + dt * k1 + dt^2 / 2 * (k2 + k3)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
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

    f(k1, uprev, p, t)
    TaylorDiff.derivative!(k2, (_y, _u) -> f(_y, _u, p, t), tmp, uprev, k1, Val(1))
    TaylorDiff.derivative!(k3, (_y, _t) -> f(_y, uprev, p, _t), tmp, t, one(t), Val(1))
    @.. u = uprev + dt * k1 + dt^2 / 2 * (k2 + k3)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end
