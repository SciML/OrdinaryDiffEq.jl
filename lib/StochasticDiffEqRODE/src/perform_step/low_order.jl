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
    f0 = integrator.f(utilde, p, t + dt, w0)
    fp = integrator.f(utilde, p, t + dt, w0 + h)
    fm = integrator.f(utilde, p, t + dt, w0 - h)
    u = uprev .+ (dt / 2) .* (ftmp .+ f0) .+ (I1 / (2 * h)) .* (fp .- fm) .+
        (I2 / (2 * adt)) .* (fp .- 2 .* f0 .+ fm)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTaylor15Cache)
    (; tmp, rtmp, rtmp0, rtmpp, rtmpm) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    adt = abs(dt)
    h = sqrt(adt)
    I1, I2 = path_integrals(W, t, dt, w0)
    integrator.f(rtmp, uprev, p, t, w0)
    @.. tmp = uprev + dt * rtmp
    integrator.f(rtmp0, tmp, p, t + dt, w0)
    integrator.f(rtmpp, tmp, p, t + dt, w0 + h)
    integrator.f(rtmpm, tmp, p, t + dt, w0 - h)
    @.. u = uprev + (dt / 2) * (rtmp + rtmp0) + (I1 / (2 * h)) * (rtmpp - rtmpm) +
        (I2 / (2 * adt)) * (rtmpp - 2 * rtmp0 + rtmpm)
end

@inline function segment_integrals(L, s0, a, b)
    a2 = a * a
    b2 = b * b
    ab = a * b
    i1 = L * (a + b) / 2
    i2 = L * (a2 + ab + b2) / 3
    i3 = L * (a + b) * (a2 + b2) / 4
    i4 = L * (a2 * a2 + ab * (a2 + b2) + a2 * b2 + b2 * b2) / 5
    m1 = s0 * i1 + L^2 * (a + 2 * b) / 6
    m2 = s0 * i2 + L^2 * (a2 + 2 * ab + 3 * b2) / 12
    return (i1, i2, i3, i4, m1, m2)
end

function path_integrals25(W, t, dt, w0)
    tend = t + dt
    tdir = dt < zero(dt) ? -one(dt) : one(dt)
    rev = length(W.t) > 1 && W.t[2] < W.t[1]
    ilo = searchsortedfirst(W.t, t, rev = rev)
    ihi = searchsortedlast(W.t, tend, rev = rev)
    acc = segment_integrals(zero(dt), zero(dt), zero(w0), zero(w0))
    tprev = t
    vprev = zero(w0)
    @inbounds for i in ilo:ihi
        ti = W.t[i]
        tdir * (ti - tprev) <= zero(dt) && continue
        vi = W.W[i] - w0
        acc = acc .+ segment_integrals(ti - tprev, tprev - t, vprev, vi)
        tprev = ti
        vprev = vi
    end
    return acc .+ segment_integrals(tend - tprev, tprev - t, vprev, W.dW)
end

@muladd function perform_step!(integrator, cache::RandomTaylor25ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    adt = abs(dt)
    δ = sqrt(adt)
    ε = dt / 2
    I1, I2, I3, I4, I10, I20 = path_integrals25(W, t, dt, w0)
    I01 = dt * I1 - I10
    I02 = dt * I2 - I20
    I11 = I1^2 / 2
    g0 = integrator.f(uprev, p, t, w0)
    gp1 = integrator.f(uprev, p, t, w0 + δ)
    gm1 = integrator.f(uprev, p, t, w0 - δ)
    gp2 = integrator.f(uprev, p, t, w0 + 2 * δ)
    gm2 = integrator.f(uprev, p, t, w0 - 2 * δ)
    gw = (8 .* (gp1 .- gm1) .- (gp2 .- gm2)) ./ (12 * δ)
    gww = (16 .* (gp1 .+ gm1) .- (gp2 .+ gm2) .- 30 .* g0) ./ (12 * adt)
    gwww = ((gp2 .- gm2) .- 2 .* (gp1 .- gm1)) ./ (2 * δ * adt)
    gwwww = ((gp2 .+ gm2) .- 4 .* (gp1 .+ gm1) .+ 6 .* g0) ./ adt^2
    u1 = uprev .+ ε .* g0
    u2 = uprev .+ dt .* g0
    f1 = integrator.f(u1, p, t + ε, w0)
    f2 = integrator.f(u2, p, t + dt, w0)
    jg = (4 .* f1 .- f2 .- 3 .* g0) ./ dt
    hgg = (g0 .- 2 .* f1 .+ f2) ./ ε^2
    jjg = (integrator.f(uprev .+ ε .* jg, p, t, w0) .- g0) ./ ε
    jgp = (4 .* integrator.f(u1, p, t + ε, w0 + δ) .- integrator.f(u2, p, t + dt, w0 + δ) .- 3 .* gp1) ./ dt
    jgm = (4 .* integrator.f(u1, p, t + ε, w0 - δ) .- integrator.f(u2, p, t + dt, w0 - δ) .- 3 .* gm1) ./ dt
    s1 = adt * δ
    up = uprev .+ s1 .* gw
    um = uprev .- s1 .* gw
    jgwp = (integrator.f(up, p, t, w0 + δ) .- integrator.f(um, p, t, w0 + δ)) ./ (2 * s1)
    jgwm = (integrator.f(up, p, t, w0 - δ) .- integrator.f(um, p, t, w0 - δ)) ./ (2 * s1)
    s2 = adt^2
    jgww = (integrator.f(uprev .+ s2 .* gww, p, t, w0) .- integrator.f(uprev .- s2 .* gww, p, t, w0)) ./ (2 * s2)
    u = uprev .+ dt .* g0 .+ I1 .* gw .+ (I2 / 2) .* gww .+ (I3 / 6) .* gwww .+ (I4 / 24) .* gwwww .+
        (dt^2 / 2) .* jg .+ (I10 / (2 * δ)) .* (jgp .- jgm) .+ (I20 / (2 * adt)) .* (jgp .- 2 .* jg .+ jgm) .+
        (I01 / 2) .* (jgwp .+ jgwm) .+ (I11 / (2 * δ)) .* (jgwp .- jgwm) .+ (I02 / 2) .* jgww .+
        (dt^3 / 6) .* (hgg .+ jjg)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTaylor25Cache)
    (; tmp, g0, gp1, gm1, gp2, gm2, gw, gww, jg, rtmp1, rtmp2, rtmp3) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    adt = abs(dt)
    δ = sqrt(adt)
    ε = dt / 2
    I1, I2, I3, I4, I10, I20 = path_integrals25(W, t, dt, w0)
    I01 = dt * I1 - I10
    I02 = dt * I2 - I20
    I11 = I1^2 / 2
    integrator.f(g0, uprev, p, t, w0)
    integrator.f(gp1, uprev, p, t, w0 + δ)
    integrator.f(gm1, uprev, p, t, w0 - δ)
    integrator.f(gp2, uprev, p, t, w0 + 2 * δ)
    integrator.f(gm2, uprev, p, t, w0 - 2 * δ)
    @.. gw = (8 * (gp1 - gm1) - (gp2 - gm2)) / (12 * δ)
    @.. gww = (16 * (gp1 + gm1) - (gp2 + gm2) - 30 * g0) / (12 * adt)
    @.. u = uprev + dt * g0 + I1 * gw + (I2 / 2) * gww +
        (I3 / 6) * ((gp2 - gm2) - 2 * (gp1 - gm1)) / (2 * δ * adt) +
        (I4 / 24) * ((gp2 + gm2) - 4 * (gp1 + gm1) + 6 * g0) / adt^2
    @.. tmp = uprev + ε * g0
    integrator.f(rtmp1, tmp, p, t + ε, w0)
    @.. tmp = uprev + dt * g0
    integrator.f(rtmp2, tmp, p, t + dt, w0)
    @.. jg = (4 * rtmp1 - rtmp2 - 3 * g0) / dt
    @.. u = u + (dt^2 / 2) * jg + (dt^3 / 6) * (g0 - 2 * rtmp1 + rtmp2) / ε^2
    @.. tmp = uprev + ε * jg
    integrator.f(rtmp1, tmp, p, t, w0)
    @.. u = u + (dt^3 / 6) * (rtmp1 - g0) / ε
    @.. tmp = uprev + ε * g0
    integrator.f(rtmp1, tmp, p, t + ε, w0 + δ)
    @.. tmp = uprev + dt * g0
    integrator.f(rtmp2, tmp, p, t + dt, w0 + δ)
    @.. rtmp3 = (4 * rtmp1 - rtmp2 - 3 * gp1) / dt
    @.. tmp = uprev + ε * g0
    integrator.f(rtmp1, tmp, p, t + ε, w0 - δ)
    @.. tmp = uprev + dt * g0
    integrator.f(rtmp2, tmp, p, t + dt, w0 - δ)
    @.. rtmp1 = (4 * rtmp1 - rtmp2 - 3 * gm1) / dt
    @.. u = u + (I10 / (2 * δ)) * (rtmp3 - rtmp1) + (I20 / (2 * adt)) * (rtmp3 - 2 * jg + rtmp1)
    s1 = adt * δ
    @.. tmp = uprev + s1 * gw
    integrator.f(rtmp1, tmp, p, t, w0 + δ)
    @.. tmp = uprev - s1 * gw
    integrator.f(rtmp2, tmp, p, t, w0 + δ)
    @.. rtmp3 = (rtmp1 - rtmp2) / (2 * s1)
    @.. tmp = uprev + s1 * gw
    integrator.f(rtmp1, tmp, p, t, w0 - δ)
    @.. tmp = uprev - s1 * gw
    integrator.f(rtmp2, tmp, p, t, w0 - δ)
    @.. rtmp1 = (rtmp1 - rtmp2) / (2 * s1)
    @.. u = u + (I01 / 2) * (rtmp3 + rtmp1) + (I11 / (2 * δ)) * (rtmp3 - rtmp1)
    s2 = adt^2
    @.. tmp = uprev + s2 * gww
    integrator.f(rtmp1, tmp, p, t, w0)
    @.. tmp = uprev - s2 * gww
    integrator.f(rtmp2, tmp, p, t, w0)
    @.. u = u + (I02 / 2) * (rtmp1 - rtmp2) / (2 * s2)
end
