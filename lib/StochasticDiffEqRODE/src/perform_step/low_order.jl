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

function add_segment!(I1, I2, L, a, b)
    @inbounds for j in eachindex(a)
        I1[j] += L * (a[j] + b[j]) / 2
        for i in eachindex(a)
            I2[i, j] += L * (2 * a[i] * a[j] + a[i] * b[j] + b[i] * a[j] + 2 * b[i] * b[j]) / 6
        end
    end
    return nothing
end

function path_integrals!(I1, I2, vprev, vcur, W, t, dt, w0)
    tend = t + dt
    tdir = dt < zero(dt) ? -one(dt) : one(dt)
    rev = length(W.t) > 1 && W.t[2] < W.t[1]
    ilo = searchsortedfirst(W.t, t, rev = rev)
    ihi = searchsortedlast(W.t, tend, rev = rev)
    fill!(I1, false)
    fill!(I2, false)
    fill!(vprev, false)
    tprev = t
    @inbounds for i in ilo:ihi
        ti = W.t[i]
        tdir * (ti - tprev) <= zero(dt) && continue
        vcur .= W.W[i] .- w0
        add_segment!(I1, I2, ti - tprev, vprev, vcur)
        tprev = ti
        vprev .= vcur
    end
    add_segment!(I1, I2, tend - tprev, vprev, W.dW)
    return nothing
end

shifted_noise(w0, s, i) = w0 .+ s .* (eachindex(w0) .== i)
function shifted_noise(w0, s, i, j)
    return w0 .+ s .* ((eachindex(w0) .== i) .| (eachindex(w0) .== j))
end

@muladd function perform_step!(integrator, cache::RandomTaylor15VectorConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    m = length(w0)
    adt = abs(dt)
    h = sqrt(adt)
    I1 = zeros(eltype(w0), m)
    I2 = zeros(eltype(w0), m, m)
    path_integrals!(I1, I2, zero(I1), zero(I1), W, t, dt, w0)
    ftmp = integrator.f(uprev, p, t, w0)
    utilde = uprev .+ dt .* ftmp
    f0 = integrator.f(utilde, p, t + dt, w0)
    u = uprev .+ (dt / 2) .* (ftmp .+ f0)
    d2 = Vector{typeof(f0)}(undef, m)
    for i in 1:m
        fp = integrator.f(utilde, p, t + dt, shifted_noise(w0, h, i))
        fm = integrator.f(utilde, p, t + dt, shifted_noise(w0, -h, i))
        d2[i] = fp .- 2 .* f0 .+ fm
        u = u .+ (I1[i] / (2 * h)) .* (fp .- fm) .+ (I2[i, i] / (2 * adt)) .* d2[i]
    end
    for j in 2:m, i in 1:(j - 1)
        fp = integrator.f(utilde, p, t + dt, shifted_noise(w0, h, i, j))
        fm = integrator.f(utilde, p, t + dt, shifted_noise(w0, -h, i, j))
        u = u .+ (I2[i, j] / (2 * adt)) .* (fp .- 2 .* f0 .+ fm .- d2[i] .- d2[j])
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTaylor15VectorCache)
    (; tmp, rtmp, rtmp0, rtmpp, rtmpm, rdiag, wtmp, vprev, vcur, I1, I2) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    w0 = W.curW
    m = length(w0)
    adt = abs(dt)
    h = sqrt(adt)
    path_integrals!(I1, I2, vprev, vcur, W, t, dt, w0)
    integrator.f(rtmp, uprev, p, t, w0)
    @.. tmp = uprev + dt * rtmp
    integrator.f(rtmp0, tmp, p, t + dt, w0)
    @.. u = uprev + (dt / 2) * (rtmp + rtmp0)
    wtmp .= w0
    for i in 1:m
        wtmp[i] = w0[i] + h
        integrator.f(rtmpp, tmp, p, t + dt, wtmp)
        wtmp[i] = w0[i] - h
        integrator.f(rtmpm, tmp, p, t + dt, wtmp)
        wtmp[i] = w0[i]
        d = rdiag[i]
        @.. d = rtmpp - 2 * rtmp0 + rtmpm
        c1 = I1[i] / (2 * h)
        c2 = I2[i, i] / (2 * adt)
        @.. u = u + c1 * (rtmpp - rtmpm) + c2 * d
    end
    for j in 2:m, i in 1:(j - 1)
        wtmp[i] = w0[i] + h
        wtmp[j] = w0[j] + h
        integrator.f(rtmpp, tmp, p, t + dt, wtmp)
        wtmp[i] = w0[i] - h
        wtmp[j] = w0[j] - h
        integrator.f(rtmpm, tmp, p, t + dt, wtmp)
        wtmp[i] = w0[i]
        wtmp[j] = w0[j]
        di = rdiag[i]
        dj = rdiag[j]
        c = I2[i, j] / (2 * adt)
        @.. u = u + c * (rtmpp - 2 * rtmp0 + rtmpm - di - dj)
    end
end
