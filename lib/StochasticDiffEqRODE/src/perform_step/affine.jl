affine_endpoint!(wtmp::Number, W) = W.curW + W.dW
function affine_endpoint!(wtmp, W)
    @.. wtmp = W.curW + W.dW
    return wtmp
end

affine_unit!(wtmp::Number, i) = one(wtmp)
function affine_unit!(wtmp, i)
    fill!(wtmp, zero(eltype(wtmp)))
    wtmp[i] = one(eltype(wtmp))
    return wtmp
end

affine_weights!(wtmp::Number, I, dt) = first(I) / dt^2
function affine_weights!(wtmp, I, dt)
    @.. wtmp = I / dt^2
    return wtmp
end

function affine_segment!(I1, I01, I2, s, ds, a, b)
    @inbounds for j in eachindex(I1)
        aj = a[j]
        bj = b[j]
        for i in eachindex(I1)
            ai = a[i]
            bi = b[i]
            I2[i, j] += ds * (
                I1[i] * (aj + bj) / 2 +
                    ds * (3 * ai * aj + bi * aj + 5 * ai * bj + 3 * bi * bj) / 24
            )
        end
        I01[j] += ds * (s * (aj + bj) / 2 + ds * (aj + 2 * bj) / 6)
    end
    @inbounds for j in eachindex(I1)
        I1[j] += ds * (a[j] + b[j]) / 2
    end
    return nothing
end

function affine_integrals!(I1, I01, I2, W, t, dt, wend)
    tend = t + dt
    rev = length(W.t) > 1 && W.t[2] < W.t[1]
    ilo = searchsortedfirst(W.t, t, rev = rev)
    ihi = searchsortedlast(W.t, tend, rev = rev)
    fill!(I1, zero(eltype(I1)))
    fill!(I01, zero(eltype(I01)))
    fill!(I2, zero(eltype(I2)))
    tprev = t
    vprev = W.curW
    @inbounds for i in ilo:ihi
        ti = W.t[i]
        vi = W.W[i]
        affine_segment!(I1, I01, I2, tprev - t, ti - tprev, vprev, vi)
        tprev = ti
        vprev = vi
    end
    affine_segment!(I1, I01, I2, tprev - t, tend - tprev, vprev, wend)
    return nothing
end

@muladd function perform_step!(integrator, cache::RandomAffineHeunConstantCache)
    (; wzero, wtmp, I1, I01, I2) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    affine_integrals!(I1, I01, I2, W, t, dt, affine_endpoint!(wtmp, W))
    ftmp = integrator.f(uprev, p, t, wzero)
    tmp = uprev .+ dt .* ftmp
    f0 = integrator.f(tmp, p, t + dt, wzero)
    fw = integrator.f(tmp, p, t + dt, affine_weights!(wtmp, I01, dt))
    u = uprev .+ (dt / 2) .* (ftmp .+ f0) .+ dt .* (fw .- f0)
    for i in eachindex(I1)
        b1 = I1[i] - (I01[i] + sum(view(I2, :, i))) / dt
        b2 = I1[i] - I01[i] / dt
        fi = integrator.f(uprev, p, t, affine_unit!(wtmp, i)) .- ftmp
        tmp = uprev .+ dt .* fi
        f0 = integrator.f(tmp, p, t, wzero)
        fw = integrator.f(tmp, p, t, affine_weights!(wtmp, view(I2, i, :), dt))
        u = u .+ b1 .* fi .+ b2 .* (f0 .- ftmp) .+ dt .* (fw .- f0)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomAffineHeunCache)
    (; tmp, rtmp1, rtmp2, rtmp3, wzero, wtmp, I1, I01, I2) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    affine_integrals!(I1, I01, I2, W, t, dt, affine_endpoint!(wtmp, W))
    integrator.f(rtmp1, uprev, p, t, wzero)
    @.. tmp = uprev + dt * rtmp1
    integrator.f(rtmp2, tmp, p, t + dt, wzero)
    integrator.f(rtmp3, tmp, p, t + dt, affine_weights!(wtmp, I01, dt))
    @.. u = uprev + (dt / 2) * (rtmp1 + rtmp2) + dt * (rtmp3 - rtmp2)
    for i in eachindex(I1)
        b1 = I1[i] - (I01[i] + sum(view(I2, :, i))) / dt
        b2 = I1[i] - I01[i] / dt
        integrator.f(rtmp2, uprev, p, t, affine_unit!(wtmp, i))
        @.. rtmp2 = rtmp2 - rtmp1
        @.. tmp = uprev + dt * rtmp2
        @.. u = u + b1 * rtmp2
        integrator.f(rtmp2, tmp, p, t, wzero)
        integrator.f(rtmp3, tmp, p, t, affine_weights!(wtmp, view(I2, i, :), dt))
        @.. u = u + b2 * (rtmp2 - rtmp1) + dt * (rtmp3 - rtmp2)
    end
end
