@muladd function perform_step!(integrator, cache::RandomEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    u = uprev .+ dt .* integrator.f(uprev, p, t, W.curW)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomEMCache)
    (; rtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(rtmp, uprev, p, t, W.curW)
    @.. u = uprev + dt * rtmp
end

@muladd function perform_step!(integrator, cache::RandomTamedEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    ftmp = integrator.f(uprev, p, t, W.curW)
    u = uprev .+ dt .* ftmp ./ (1 .+ dt .* norm(ftmp))
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTamedEMCache)
    (; rtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(rtmp, uprev, p, t, W.curW)
    tamed = 1 + dt * norm(rtmp)
    @.. u = uprev + dt * rtmp / tamed
end

@muladd function perform_step!(integrator, cache::RandomHeunConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    ftmp = integrator.f(uprev, p, t, W.curW)
    tmp = @.. uprev + dt * ftmp
    wtmp = @.. W.curW + W.dW
    u = uprev .+ (dt / 2) .* (ftmp .+ integrator.f(tmp, p, t + dt, wtmp))
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomHeunCache)
    (; tmp, rtmp1, rtmp2, wtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(rtmp1, uprev, p, t, W.curW)
    @.. tmp = uprev + dt * rtmp1
    if W.dW isa Number
        wtmp = W.curW + W.dW
    else
        @.. wtmp = W.curW + W.dW
    end
    integrator.f(rtmp2, tmp, p, t + dt, wtmp)
    @.. u = uprev + (dt / 2) * (rtmp1 + rtmp2)
end

function path_integrals(W, t, dt, w0)
    tend = t + dt
    tdir = dt < zero(dt) ? -one(dt) : one(dt)
    rev = length(W.t) > 1 && W.t[2] < W.t[1]
    ilo = searchsortedfirst(W.t, t, rev = rev)
    ihi = searchsortedlast(W.t, tend, rev = rev)
    I1 = zero(w0) * dt
    I2 = zero(w0) * zero(w0) * dt
    tprev = t
    vprev = zero(w0)
    @inbounds for i in ilo:ihi
        ti = W.t[i]
        tdir * (ti - tprev) <= zero(dt) && continue
        vi = W.W[i] - w0
        I1 += (ti - tprev) * (vprev + vi) / 2
        I2 += (ti - tprev) * (vprev^2 + vprev * vi + vi^2) / 3
        tprev = ti
        vprev = vi
    end
    vend = W.dW
    return I1 + (tend - tprev) * (vprev + vend) / 2,
        I2 + (tend - tprev) * (vprev^2 + vprev * vend + vend^2) / 3
end

@muladd function perform_step!(integrator, cache::RandomTaylor15ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    adt = abs(dt)
    h = sqrt(adt)
    I1, I2 = path_integrals(W, t, dt, w0)
    ftmp = integrator.f(uprev, p, t, w0)
    utilde = uprev .+ dt .* ftmp
    fp = integrator.f(utilde, p, t + dt, w0 + h)
    fm = integrator.f(utilde, p, t + dt, w0 - h)
    u = uprev .+ dt .* ftmp .+ (I1 / (2 * h)) .* (fp .- fm) .+
        (I2 / (2 * adt)) .* (fp .- 2 .* ftmp .+ fm)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTaylor15Cache)
    (; tmp, rtmp, rtmpp, rtmpm) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    adt = abs(dt)
    h = sqrt(adt)
    I1, I2 = path_integrals(W, t, dt, w0)
    integrator.f(rtmp, uprev, p, t, w0)
    @.. tmp = uprev + dt * rtmp
    integrator.f(rtmpp, tmp, p, t + dt, w0 + h)
    integrator.f(rtmpm, tmp, p, t + dt, w0 - h)
    @.. u = uprev + dt * rtmp + (I1 / (2 * h)) * (rtmpp - rtmpm) +
        (I2 / (2 * adt)) * (rtmpp - 2 * rtmp + rtmpm)
end
