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

function path_integral(W, t, dt, w0)
    tend = t + dt
    integral = zero(w0) * dt
    tprev = t
    vprev = zero(w0)
    @inbounds for i in searchsortedfirst(W.t, t):searchsortedlast(W.t, tend)
        ti = W.t[i]
        ti <= tprev && continue
        vi = W.W[i] - w0
        integral += (ti - tprev) * (vprev + vi) / 2
        tprev = ti
        vprev = vi
    end
    return integral + (tend - tprev) * (vprev + W.dW) / 2
end

@muladd function perform_step!(integrator, cache::RandomTaylor15ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    sqdt = sqrt(dt)
    I10 = path_integral(W, t, dt, w0)
    ftmp = integrator.f(uprev, p, t, w0)
    utilde = uprev .+ dt .* ftmp
    fp = integrator.f(utilde, p, t + dt, w0 + sqdt)
    fm = integrator.f(utilde, p, t + dt, w0 - sqdt)
    u = uprev .+ dt .* ftmp .+ (I10 / (2 * sqdt)) .* (fp .- fm) .+
        (dt / 4) .* (fp .- 2 .* ftmp .+ fm)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTaylor15Cache)
    (; tmp, rtmp, rtmpp, rtmpm) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    sqdt = sqrt(dt)
    I10 = path_integral(W, t, dt, w0)
    integrator.f(rtmp, uprev, p, t, w0)
    @.. tmp = uprev + dt * rtmp
    integrator.f(rtmpp, tmp, p, t + dt, w0 + sqdt)
    integrator.f(rtmpm, tmp, p, t + dt, w0 - sqdt)
    @.. u = uprev + dt * rtmp + (I10 / (2 * sqdt)) * (rtmpp - rtmpm) +
        (dt / 4) * (rtmpp - 2 * rtmp + rtmpm)
end
