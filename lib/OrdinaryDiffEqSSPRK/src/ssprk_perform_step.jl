function initialize!(integrator, cache::SSPRK22ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK22ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator

    # u1 -> stored as u
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + dt * integrator.fsalfirst
    k = f(u, p, t + dt)
    # u
    u = (uprev + u + dt * k) / 2

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK22Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK22Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, stage_limiter!, step_limiter!, thread) = cache

    # u1 -> stored as u
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+dt*fsalfirst
    stage_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    # u
    @.. broadcast=false thread=thread u=(uprev+u+dt*k)/2
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
end

function initialize!(integrator, cache::KYKSSPRK42Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator, cache::KYKSSPRK42Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, tmp, fsalfirst, stage_limiter!, step_limiter!, thread) = cache
    (; α20, α21, α30, α32, α40, α43, β10, β21, β30, β32, β40, β43, c1, c2, c3) = cache.tab

    δ = fsalfirst
    # u1 -> stored as u
    @.. broadcast=false thread=thread u=uprev+dt*β10*δ
    stage_limiter!(u, integrator, p, t + c1 * dt)
    f(k, u, p, t + c1 * dt)
    # u2
    @.. broadcast=false thread=thread tmp=α20*uprev+α21*u+dt*β21*k
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k, tmp, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread tmp=α30*uprev+α32*tmp+dt*β30*δ+
                                          dt*β32*k
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k, tmp, p, t + c3 * dt)
    # u
    @.. broadcast=false thread=thread u=α40*uprev+α43*tmp+dt*β40*δ+
                                        dt*β43*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
end

function initialize!(integrator, cache::KYKSSPRK42ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::KYKSSPRK42ConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; α20, α21, α30, α32, α40, α43, β10, β21, β30, β32, β40, β43, c1, c2, c3) = cache

    #u1
    δ = integrator.fsalfirst
    u = uprev + dt * β10 * δ
    k = f(u, p, t + c1 * dt)
    #u2
    tmp = α20 * uprev + α21 * u + dt * β21 * k
    k = f(tmp, p, t + c2 * dt)
    #u3
    tmp = α30 * uprev + α32 * tmp + dt * β30 * δ + dt * β32 * k
    k = f(tmp, p, t + c3 * dt)
    #u
    u = α40 * uprev + α43 * tmp + dt * β40 * δ + dt * β43 * k

    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    integrator.k[1] = integrator.fsalfirst
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK33ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK33ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + dt * integrator.fsalfirst
    k = f(u, p, t + dt)
    # u2
    u = (3 * uprev + u + dt * k) / 4
    k = f(u, p, t + dt / 2)
    # u
    u = (uprev + 2 * u + 2 * dt * k) / 3

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK33Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK33Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, stage_limiter!, step_limiter!, thread) = cache

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+dt*fsalfirst
    stage_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    # u2
    @.. broadcast=false thread=thread u=(3*uprev+u+dt*k)/4
    stage_limiter!(u, integrator, p, t + dt / 2)
    f(k, u, p, t + dt / 2)
    # u
    @.. broadcast=false thread=thread u=(uprev+2*u+2*dt*k)/3
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
end

function initialize!(integrator, cache::SSPRK53ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α30, α32, α40, α43, α52, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    tmp = uprev + β10 * dt * integrator.fsalfirst
    k = f(tmp, p, t + c1 * dt)
    # u2 -> stored as u
    u = tmp + β21 * dt * k
    k = f(u, p, t + c2 * dt)
    # u3
    tmp = α30 * uprev + α32 * u + β32 * dt * k
    k = f(tmp, p, t + c3 * dt)
    # u4
    tmp = α40 * uprev + α43 * tmp + β43 * dt * k
    k = f(tmp, p, t + c4 * dt)
    # u
    u = α52 * u + α54 * tmp + β54 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK53Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; α30, α32, α40, α43, α52, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache.tab

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread tmp=uprev+β10*dt*fsalfirst
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k, tmp, p, t + c1 * dt)
    # u2 -> stored as u
    @.. broadcast=false thread=thread u=tmp+β21*dt*k
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(k, u, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread tmp=α30*uprev+α32*u+β32*dt*k
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k, tmp, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread tmp=α40*uprev+α43*tmp+β43*dt*k
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k, tmp, p, t + c4 * dt)
    # u
    @.. broadcast=false thread=thread u=α52*u+α54*tmp+β54*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
end

function initialize!(integrator, cache::SSPRK53_2N1ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53_2N1ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α40, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache
    #stores in u for all intermediate stages
    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + β10 * dt * integrator.fsalfirst
    k = f(u, p, t + c1 * dt)
    # u2
    u = u + β21 * dt * k
    k = f(u, p, t + c2 * dt)
    # u3
    u = u + β32 * dt * k
    k = f(u, p, t + c3 * dt)
    # u4
    u = α40 * uprev + α43 * u + β43 * dt * k
    k = f(u, p, t + c4 * dt)
    # u
    u = u + β54 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK53_2N1Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53_2N1Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, stage_limiter!, step_limiter!, thread) = cache
    (; α40, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache.tab

    #stores in u for all intermediate stages
    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+β10*dt*fsalfirst
    stage_limiter!(u, integrator, p, t + c1 * dt)
    f(k, u, p, t + c1 * dt)
    # u2
    @.. broadcast=false thread=thread u=u+β21*dt*k
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(k, u, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread u=u+β32*dt*k
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(k, u, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread u=α40*uprev+α43*u+β43*dt*k
    stage_limiter!(u, integrator, p, t + c4 * dt)
    f(k, u, p, t + c4 * dt)
    # u
    @.. broadcast=false thread=thread u=u+β54*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
end

function initialize!(integrator, cache::SSPRK53_2N2ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53_2N2ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α30, α32, α50, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache
    #stores in u for all intermediate stages
    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + β10 * dt * integrator.fsalfirst
    k = f(u, p, t + c1 * dt)
    # u2 -> stored as u
    u = u + β21 * dt * k
    k = f(u, p, t + c2 * dt)
    # u3
    u = α30 * uprev + α32 * u + β32 * dt * k
    k = f(u, p, t + c3 * dt)
    # u4
    u = u + β43 * dt * k
    k = f(u, p, t + c4 * dt)
    # u
    u = α50 * uprev + α54 * u + β54 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK53_2N2Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53_2N2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, stage_limiter!, step_limiter!, thread) = cache
    (; α30, α32, α50, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache.tab

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+β10*dt*fsalfirst
    stage_limiter!(u, integrator, p, t + c1 * dt)
    f(k, u, p, t + c1 * dt)
    # u2 -> stored as u
    @.. broadcast=false thread=thread u=u+β21*dt*k
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(k, u, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread u=α30*uprev+α32*u+β32*dt*k
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(k, u, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread u=u+β43*dt*k
    stage_limiter!(u, integrator, p, t + c4 * dt)
    f(k, u, p, t + c4 * dt)
    # u
    @.. broadcast=false thread=thread u=α50*uprev+α54*u+β54*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
end

function initialize!(integrator, cache::SSPRK53_HConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53_HConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α30, α32, α40, α41, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache
    #stores in u for all intermediate stages
    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    tmp = uprev + β10 * dt * integrator.fsalfirst
    k = f(tmp, p, t + c1 * dt)
    # u2
    u = tmp + β21 * dt * k
    k = f(u, p, t + c2 * dt)
    # u3
    u = α30 * uprev + α32 * u + β32 * dt * k
    k = f(u, p, t + c3 * dt)
    # u4
    u = α40 * uprev + α41 * tmp + α43 * u + β43 * dt * k
    k = f(u, p, t + c4 * dt)
    # u
    u = u + β54 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK53_HCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK53_HCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; α30, α32, α40, α41, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4) = cache.tab
    #stores in u for all intermediate stages
    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread tmp=uprev+β10*dt*fsalfirst
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k, tmp, p, t + c1 * dt)
    # u2
    @.. broadcast=false thread=thread u=tmp+β21*dt*k
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(k, u, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread u=α30*uprev+α32*u+β32*dt*k
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(k, u, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread u=α40*uprev+α41*tmp+α43*u+β43*dt*k
    stage_limiter!(u, integrator, p, t + c4 * dt)
    f(k, u, p, t + c4 * dt)
    # u
    @.. broadcast=false thread=thread u=u+β54*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
end

function initialize!(integrator, cache::SSPRK63ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK63ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4, c5) = cache

    # u1 -> stored as u
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + β10 * dt * integrator.fsalfirst
    k = f(u, p, t + c1 * dt)
    # u2
    u₂ = u + β21 * dt * k
    k = f(u₂, p, t + c2 * dt)
    # u3
    tmp = u₂ + β32 * dt * k
    k = f(tmp, p, t + c3 * dt)
    # u4
    tmp = α40 * uprev + α41 * u + α43 * tmp + β43 * dt * k
    k = f(tmp, p, t + c4 * dt)
    # u5
    tmp = tmp + β54 * dt * k
    k = f(tmp, p, t + c5 * dt)
    # u
    u = α62 * u₂ + α65 * tmp + β65 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK63Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK63Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, u₂, stage_limiter!, step_limiter!, thread) = cache
    (; α40, α41, α43, α62, α65, β10, β21, β32, β43,
        β54, β65, c1, c2, c3, c4, c5) = cache.tab

    # u1 -> stored as u
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+β10*dt*fsalfirst
    stage_limiter!(u, integrator, p, t + c1 * dt)
    f(k, u, p, t + c1 * dt)
    # u2
    @.. broadcast=false thread=thread u₂=u+β21*dt*k
    stage_limiter!(u₂, integrator, p, t + c2 * dt)
    f(k, u₂, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread tmp=u₂+β32*dt*k
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k, tmp, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread tmp=α40*uprev+α41*u+α43*tmp+β43*dt*k
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k, tmp, p, t + c4 * dt)
    # u5
    @.. broadcast=false thread=thread tmp=tmp+β54*dt*k
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k, tmp, p, t + c5 * dt)
    # u
    @.. broadcast=false thread=thread u=α62*u₂+α65*tmp+β65*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
end

function initialize!(integrator, cache::SSPRK73ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK73ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α40, α43, α50, α51, α54, α73, α76, β10, β21, β32,
        β43, β54, β65, β76, c1, c2, c3, c4, c5, c6) = cache

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u₁ = uprev + β10 * dt * integrator.fsalfirst
    k = f(u₁, p, t + c1 * dt)
    # u2
    tmp = u₁ + β21 * dt * k
    k = f(tmp, p, t + c2 * dt)
    # u3 -> stored as u
    u = tmp + β32 * dt * k
    k = f(u, p, t + c3 * dt)
    # u4
    tmp = α40 * uprev + α43 * u + β43 * dt * k
    k = f(tmp, p, t + c4 * dt)
    # u5
    tmp = α50 * uprev + α51 * u₁ + α54 * tmp + β54 * dt * k
    k = f(tmp, p, t + c5 * dt)
    # u6
    tmp = tmp + β65 * dt * k
    k = f(tmp, p, t + c6 * dt)
    # u
    u = α73 * u + α76 * tmp + β76 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK73Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK73Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, u₁, stage_limiter!, step_limiter!, thread) = cache
    (; α40, α43, α50, α51, α54, α73, α76, β10, β21, β32, β43,
        β54, β65, β76, c1, c2, c3, c4, c5, c6) = cache.tab

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u₁=uprev+β10*dt*fsalfirst
    stage_limiter!(u₁, integrator, p, t + c1 * dt)
    f(k, u₁, p, t + c1 * dt)
    # u2
    @.. broadcast=false thread=thread tmp=u₁+β21*dt*k
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k, tmp, p, t + c2 * dt)
    # u3 -> stored as u
    @.. broadcast=false thread=thread u=tmp+β32*dt*k
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(k, u, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread tmp=α40*uprev+α43*u+β43*dt*k
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k, tmp, p, t + c4 * dt)
    # u5
    @.. broadcast=false thread=thread tmp=α50*uprev+α51*u₁+α54*tmp+β54*dt*k
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k, tmp, p, t + c5 * dt)
    # u6
    @.. broadcast=false thread=thread tmp=tmp+β65*dt*k
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k, tmp, p, t + c6 * dt)
    # u
    @.. broadcast=false thread=thread u=α73*u+α76*tmp+β76*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
end

function initialize!(integrator, cache::SSPRK83ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK83ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α50, α51, α54, α61, α65, α72, α73, α76, β10, β21, β32, β43,
        β54, β65, β76, β87, c1, c2, c3, c4, c5, c6, c7) = cache

    # u1 -> save as u
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + β10 * dt * integrator.fsalfirst
    k = f(u, p, t + c1 * dt)
    # u2
    u₂ = u + β21 * dt * k
    k = f(u₂, p, t + c2 * dt)
    # u3
    u₃ = u₂ + β32 * dt * k
    k = f(u₃, p, t + c3 * dt)
    # u4
    tmp = u₃ + β43 * dt * k
    k = f(tmp, p, t + c4 * dt)
    # u5
    tmp = α50 * uprev + α51 * u + α54 * tmp + β54 * dt * k
    k = f(tmp, p, t + c5 * dt)
    # u6
    tmp = α61 * u + α65 * tmp + β65 * dt * k
    k = f(tmp, p, t + c6 * dt)
    # u7
    tmp = α72 * u₂ + α73 * u₃ + α76 * tmp + β76 * dt * k
    k = f(tmp, p, t + c7 * dt)
    # u
    u = tmp + β87 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK83Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK83Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, u₂, u₃, stage_limiter!, step_limiter!, thread) = cache
    (; α50, α51, α54, α61, α65, α72, α73, α76, β10, β21, β32, β43,
        β54, β65, β76, β87, c1, c2, c3, c4, c5, c6, c7) = cache.tab

    # u1 -> save as u
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+β10*dt*fsalfirst
    stage_limiter!(u, integrator, p, t + c1 * dt)
    f(k, u, p, t + c1 * dt)
    # u2
    @.. broadcast=false thread=thread u₂=u+β21*dt*k
    stage_limiter!(u₂, integrator, p, t + c2 * dt)
    f(k, u₂, p, t + c2 * dt)
    # u3
    @.. broadcast=false thread=thread u₃=u₂+β32*dt*k
    stage_limiter!(u₃, integrator, p, t + c3 * dt)
    f(k, u₃, p, t + c3 * dt)
    # u4
    @.. broadcast=false thread=thread tmp=u₃+β43*dt*k
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k, tmp, p, t + c4 * dt)
    # u5
    @.. broadcast=false thread=thread tmp=α50*uprev+α51*u+α54*tmp+β54*dt*k
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k, tmp, p, t + c5 * dt)
    # u6
    @.. broadcast=false thread=thread tmp=α61*u+α65*tmp+β65*dt*k
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k, tmp, p, t + c6 * dt)
    # u7
    @.. broadcast=false thread=thread tmp=α72*u₂+α73*u₃+α76*tmp+β76*dt*k
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k, tmp, p, t + c7 * dt)
    # u
    @.. broadcast=false thread=thread u=tmp+β87*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)
end

function initialize!(integrator, cache::SSPRK43ConstantCache)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK43ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; one_third_u, two_thirds_u, half_u, half_t) = cache
    dt_2 = half_t * dt

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + dt_2 * integrator.fsalfirst
    k = f(u, p, t + dt_2)
    # u2
    u = u + dt_2 * k
    k = f(u, p, t + dt)
    u = u + dt_2 * k
    utilde = u  # Initialize for JET
    if integrator.opts.adaptive
        utilde = one_third_u * uprev + two_thirds_u * u # corresponds to bhat = (1/3, 1/3, 1/3, 0)
    end
    # u3
    u = two_thirds_u * uprev + one_third_u * u
    k = f(u, p, t + dt_2)
    # u
    u = u + dt_2 * k # corresponds to b = (1/6, 1/6, 1/6, 1/2)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    if integrator.opts.adaptive
        utilde = half_u * (utilde - u) # corresponds to bhat = (1/4, 1/4, 1/4, 1/4)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK43Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK43Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; one_third_u, two_thirds_u, half_u, half_t) = cache.tab
    dt_2 = half_t * dt

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+dt_2*fsalfirst
    stage_limiter!(u, integrator, p, t + dt_2)
    f(k, u, p, t + dt_2)
    # u2
    @.. broadcast=false thread=thread u=u+dt_2*k
    stage_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    #
    @.. broadcast=false thread=thread u=u+dt_2*k
    stage_limiter!(u, integrator, p, t + dt + dt_2)
    if integrator.opts.adaptive
        @.. broadcast=false utilde=one_third_u*uprev+two_thirds_u*u # corresponds to bhat = (1/3, 1/3, 1/3, 0)
    end
    # u3
    @.. broadcast=false thread=thread u=two_thirds_u*uprev+one_third_u*u
    f(k, u, p, t + dt_2)
    #
    @.. broadcast=false thread=thread u=u+dt_2*k # corresponds to b = (1/6, 1/6, 1/6, 1/2)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=half_u*(utilde-u) # corresponds to bhat = (1/4, 1/4, 1/4, 1/4)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
end

function initialize!(integrator, cache::SSPRK432ConstantCache)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK432ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    dt_2 = dt / 2

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + dt_2 * integrator.fsalfirst
    k = f(u, p, t + dt_2)
    # u2
    u = u + dt_2 * k
    k = f(u, p, t + dt)
    u = u + dt_2 * k
    utilde = u  # Initialize for JET
    if integrator.opts.adaptive
        utilde = (uprev + 2 * u) / 3
    end
    # u3
    u = (2 * uprev + u) / 3
    k = f(u, p, t + dt_2)
    # u
    u = u + dt_2 * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    if integrator.opts.adaptive
        utilde = utilde - u
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK432Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK432Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    dt_2 = dt / 2

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+dt_2*fsalfirst
    stage_limiter!(u, integrator, p, t + dt_2)
    f(k, u, p, t + dt_2)
    # u2
    @.. broadcast=false thread=thread u=u+dt_2*k
    stage_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    #
    @.. broadcast=false thread=thread u=u+dt_2*k
    stage_limiter!(u, integrator, p, t + dt + dt_2)
    if integrator.opts.adaptive
        @.. broadcast=false utilde=(uprev+2*u)/3
    end
    # u3
    @.. broadcast=false thread=thread u=(2*uprev+u)/3
    f(k, u, p, t + dt_2)
    #
    @.. broadcast=false thread=thread u=u+dt_2*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=utilde-u
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
end

function initialize!(integrator, cache::SSPRKMSVS32ConstantCache)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRKMSVS32ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; u_1, u_2, dts, dtf, μ, v_n) = cache

    if integrator.iter == 1
        cache.dts[1] = dt
        cache.dts[2] = dt
        cache.dtf[1] = dt
    end
    accpt = true
    dt = dts[1]

    if cache.step < 3 #starting Procedure
        k = f(u, p, t + dt)
        u = uprev + dt * k
        k = f(u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
        u = (uprev + u + dt * k) / 2
        if cache.step == 1
            u_2 = uprev
        else
            u_1 = uprev
        end
        if integrator.opts.adaptive
            v_n = dt / dts[2] * 0.5
            cache.dtf[2] = dtf[1]
            cache.dtf[1] = dt / v_n * 0.5
            if v_n > 0.5
                cache.step -= 1
                accpt = false
            end
            cache.dts[3] = dts[2]
            cache.dts[2] = dt
            dt = 0.9 * dtf[1]
            μ = min(dtf[1], dtf[2])
        end
    else
        if integrator.opts.adaptive
            Ω = (dts[2] + dts[3]) / dt
        else
            Ω = 2
        end
        u = (Ω * Ω - 1) / (Ω * Ω) * (uprev + Ω / (Ω - 1) * dt * integrator.fsalfirst) +
            1 / (Ω * Ω) * u_2
        u_2 = u_1
        u_1 = uprev
        if integrator.opts.adaptive
            v_n = (dts[2] + dts[3] - dt) / (dts[2] + dts[3]) * 0.5
            dt = (dts[2] + dts[3]) / (dts[2] + dts[3] + μ) * μ
            cache.dtf[2] = dtf[1]
            dtf[1] = dt / v_n * 0.5
            cache.dts[3] = dts[2]
            cache.dts[2] = dt
            μ = min(dtf[1], dtf[2])
        end
    end
    if accpt == true
        integrator.fsallast = f(u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        integrator.k[1] = integrator.fsalfirst
        integrator.u = u
    else
        integrator.fsallast = f(uprev, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        integrator.k[1] = integrator.fsalfirst
        integrator.u = uprev
    end
    cache.dts[1] = dt
    cache.step += 1
    cache.u_1 = u_1
    cache.u_2 = u_2
    cache.μ = μ
end

function initialize!(integrator, cache::SSPRKMSVS32Cache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
    integrator.fsallast = cache.k
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::SSPRKMSVS32Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, u_1, u_2, stage_limiter!, step_limiter!, thread) = cache

    if cache.step < 3
        @.. broadcast=false thread=thread u=uprev+dt*fsalfirst
        stage_limiter!(u, integrator, p, t + dt)
        f(k, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false thread=thread u=(uprev+u+dt*k)/2
        stage_limiter!(u, integrator, p, t + dt)
        step_limiter!(u, integrator, p, t + dt)

        if cache.step == 1
            cache.u_2 .= uprev
        else
            cache.u_1 .= uprev
        end
    else
        Ω = 2
        @.. broadcast=false thread=thread u=((Ω*Ω-1)/(Ω*Ω))*
        (uprev+(Ω/(Ω-1))*dt*fsalfirst)+
        (1/(Ω*Ω))*cache.u_2
        cache.u_2 .= u_1
        cache.u_1 .= uprev
        stage_limiter!(u, integrator, p, t + dt)
        step_limiter!(u, integrator, p, t + dt)
    end
    cache.step += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(k, u, p, t + dt)
end

function initialize!(integrator, cache::SSPRKMSVS43ConstantCache)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRKMSVS43ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; u_1, u_2, u_3, k1, k2, k3) = cache

    if cache.step < 4
        u = uprev + dt * integrator.fsalfirst
        k = f(u, p, t + dt)
        u = (uprev + u + dt * k) / 2
        if cache.step == 1
            u_3 = uprev
            cache.k3 = f(u_3, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        if cache.step == 2
            u_2 = uprev
            cache.k2 = f(u_2, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        if cache.step == 3
            u_1 = uprev
            cache.k1 = f(u_1, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        # u
    else
        u = (16 / 27) * (uprev + 3 * dt * integrator.fsalfirst) +
            (11 / 27) * (u_3 + (12 / 11) * dt * k3)
        cache.k3 = k2
        cache.k2 = k1
        cache.k1 = integrator.fsalfirst
        u_3 = u_2
        u_2 = u_1
        u_1 = uprev
    end
    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.u = u
    cache.step += 1
    cache.u_1 = u_1
    cache.u_2 = u_2
    cache.u_3 = u_3
end

function initialize!(integrator, cache::SSPRKMSVS43Cache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
    integrator.fsallast = cache.k
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::SSPRKMSVS43Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, u_1, u_2, u_3, stage_limiter!,
        step_limiter!, thread, k1, k2, k3) = cache

    if cache.step < 4
        @.. broadcast=false thread=thread u=uprev+dt*fsalfirst
        stage_limiter!(u, integrator, p, t + dt)
        f(k, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false thread=thread u=(uprev+u+dt*k)/2
        stage_limiter!(u, integrator, p, t + dt)
        step_limiter!(u, integrator, p, t + dt)
        if cache.step == 1
            cache.u_3 .= uprev
            f(k3, u_3, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        if cache.step == 2
            cache.u_2 .= uprev
            f(k2, u_2, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        if cache.step == 3
            cache.u_1 .= uprev
            f(k1, u_1, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        # u
    else
        @.. broadcast=false thread=thread u=(16/27)*(uprev+3*dt*fsalfirst)+
        (11/27)*(u_3+(12/11)*dt*k3)
        stage_limiter!(u, integrator, p, t + dt)
        step_limiter!(u, integrator, p, t + dt)
        cache.k3 .= k2
        cache.k2 .= k1
        cache.k1 .= fsalfirst
        cache.u_3 .= u_2
        cache.u_2 .= u_1
        cache.u_1 .= uprev
    end
    cache.step += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f(k, u, p, t + dt)
end

function initialize!(integrator, cache::SSPRK932ConstantCache)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK932ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    dt_6 = dt / 6
    dt_3 = dt / 3
    dt_2 = dt / 2

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u = uprev + dt_6 * integrator.fsalfirst
    k = f(u, p, t + dt_6)
    # u2
    u = u + dt_6 * k
    k = f(u, p, t + dt_3)
    # u3
    u = u + dt_6 * k
    k = f(u, p, t + dt_2)
    # u4
    u = u + dt_6 * k
    k = f(u, p, t + 2 * dt_3)
    # u5
    u = u + dt_6 * k
    k = f(u, p, t + 5 * dt_6)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    # u6
    u = u + dt_6 * k
    utilde = u  # Initialize for JET
    if integrator.opts.adaptive
        k = f(u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        utilde = (uprev + 6 * u + 6 * dt * k) / 7
    end
    # u6*
    u = (3 * uprev + dt_2 * integrator.fsalfirst + 2 * u) / 5
    k = f(u, p, t + dt_2)
    # u7*
    u = u + dt_6 * k
    k = f(u, p, t + 2 * dt_3)
    # u8*
    u = u + dt_6 * k
    k = f(u, p, t + 5 * dt_6)
    # u
    u = u + dt_6 * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    if integrator.opts.adaptive
        utilde = utilde - u
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK932Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK932Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    dt_6 = dt / 6
    dt_3 = dt / 3
    dt_2 = dt / 2

    # u1
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u=uprev+dt_6*fsalfirst
    stage_limiter!(u, integrator, p, t + dt_6)
    f(k, u, p, t + dt_6)
    # u2
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + dt_3)
    f(k, u, p, t + dt_3)
    # u3
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + dt_2)
    f(k, u, p, t + dt_2)
    # u4
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + 2 * dt_3)
    f(k, u, p, t + 2 * dt_3)
    # u5
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + 5 * dt_6)
    f(k, u, p, t + 5 * dt_6)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    # u6
    @.. broadcast=false thread=thread u=u+dt_6*k
    if integrator.opts.adaptive
        stage_limiter!(u, integrator, p, t + dt)
        f(k, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false thread=thread utilde=(uprev+6*u+6*dt*k)/7
    end
    # u6*
    @.. broadcast=false thread=thread u=(3*uprev+dt_2*integrator.fsalfirst+2*u)/
    5
    stage_limiter!(u, integrator, p, t + dt_6)
    f(k, u, p, t + dt_2)
    # u7*
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + 2 * dt_3)
    f(k, u, p, t + 2 * dt_3)
    # u8*
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + 5 * dt_6)
    f(k, u, p, t + 5 * dt_6)
    # u9*
    @.. broadcast=false thread=thread u=u+dt_6*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)

    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=utilde-u
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::SSPRK54ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK54ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; β10, α20, α21, β21, α30, α32, β32, α40, α43, β43,
        α52, α53, β53, α54, β54, c1, c2, c3, c4) = cache

    # u₁
    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    u₂ = uprev + β10 * dt * integrator.fsalfirst
    k = f(u₂, p, t + c1 * dt)
    # u₂
    u₂ = α20 * uprev + α21 * u₂ + β21 * dt * k
    k = f(u₂, p, t + c2 * dt)
    # u₃
    u₃ = α30 * uprev + α32 * u₂ + β32 * dt * k
    k₃ = f(u₃, p, t + c3 * dt)
    # u₄ -> stored as tmp
    tmp = α40 * uprev + α43 * u₃ + β43 * dt * k₃
    k = f(tmp, p, t + c4 * dt)
    # u
    u = α52 * u₂ + α53 * u₃ + β53 * dt * k₃ + α54 * tmp + β54 * dt * k

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK54Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK54Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, k₃, u₂, u₃, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; β10, α20, α21, β21, α30, α32, β32, α40, α43, β43,
        α52, α53, β53, α54, β54, c1, c2, c3, c4) = cache.tab

    # u₁
    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread u₂=uprev+β10*dt*fsalfirst
    stage_limiter!(u₂, integrator, p, t + c1 * dt)
    f(k, u₂, p, t + c1 * dt)
    # u₂
    @.. broadcast=false thread=thread u₂=α20*uprev+α21*u₂+β21*dt*k
    stage_limiter!(u₂, integrator, p, t + c2 * dt)
    f(k, u₂, p, t + c2 * dt)
    # u₃
    @.. broadcast=false thread=thread u₃=α30*uprev+α32*u₂+β32*dt*k
    stage_limiter!(u₃, integrator, p, t + c3 * dt)
    f(k₃, u₃, p, t + c3 * dt)
    # u₄ -> stored as tmp
    @.. broadcast=false thread=thread tmp=α40*uprev+α43*u₃+β43*dt*k₃
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k, tmp, p, t + c4 * dt)
    # u
    @.. broadcast=false thread=thread u=α52*u₂+α53*u₃+β53*dt*k₃+α54*tmp+
                                        β54*dt*k
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
end

function initialize!(integrator, cache::SSPRK104ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK104ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    dt_6 = dt / 6
    dt_3 = dt / 3
    dt_2 = dt / 2

    integrator.fsalfirst = f(uprev, p, t)
    integrator.k[1] = integrator.fsalfirst
    tmp = uprev + dt_6 * integrator.fsalfirst # u₁
    k = f(tmp, p, t + dt_6)
    tmp = tmp + dt_6 * k # u₂
    k = f(tmp, p, t + dt_3)
    tmp = tmp + dt_6 * k # u₃
    k = f(tmp, p, t + dt_2)
    u = tmp + dt_6 * k # u₄
    k₄ = f(u, p, t + 2 * dt_3)
    tmp = (3 * uprev + 2 * u + 2 * dt_6 * k₄) / 5 # u₅
    k = f(tmp, p, t + dt_3)
    tmp = tmp + dt_6 * k # u₆
    k = f(tmp, p, t + dt_2)
    tmp = tmp + dt_6 * k # u₇
    k = f(tmp, p, t + 2 * dt_3)
    tmp = tmp + dt_6 * k # u₈
    k = f(tmp, p, t + 5 * dt_6)
    tmp = tmp + dt_6 * k # u₉
    k = f(tmp, p, t + dt)
    u = (uprev + 9 * (u + dt_6 * k₄) + 15 * (tmp + dt_6 * k)) / 25

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 10)
    integrator.u = u
end

function initialize!(integrator, cache::SSPRK104Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::SSPRK104Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, k₄, tmp, stage_limiter!, step_limiter!, thread) = cache
    dt_6 = dt / 6
    dt_3 = dt / 3
    dt_2 = dt / 2

    f(fsalfirst, uprev, p, t)
    @.. broadcast=false thread=thread tmp=uprev+dt_6*fsalfirst
    stage_limiter!(tmp, integrator, p, t + dt_6)
    f(k, tmp, p, t + dt_6)
    @.. broadcast=false tmp=tmp+dt_6*k
    stage_limiter!(tmp, integrator, p, t + dt_3)
    f(k, tmp, p, t + dt_3)
    @.. broadcast=false thread=thread tmp=tmp+dt_6*k
    stage_limiter!(tmp, integrator, p, t + dt_2)
    f(k, tmp, p, t + dt_2)
    @.. broadcast=false thread=thread u=tmp+dt_6*k
    stage_limiter!(u, integrator, p, t + 2 * dt_3)
    f(k₄, u, p, t + 2 * dt_3)
    @.. broadcast=false thread=thread tmp=(3*uprev+2*u+2*dt_6*k₄)/5
    stage_limiter!(tmp, integrator, p, t + dt_3)
    f(k, tmp, p, t + dt_3)
    @.. broadcast=false thread=thread tmp=tmp+dt_6*k
    stage_limiter!(tmp, integrator, p, t + dt_2)
    f(k, tmp, p, t + dt_2)
    @.. broadcast=false thread=thread tmp=tmp+dt_6*k
    stage_limiter!(tmp, integrator, p, t + 2 * dt_3)
    f(k, tmp, p, t + 2 * dt_3)
    @.. broadcast=false thread=thread tmp=tmp+dt_6*k
    stage_limiter!(tmp, integrator, p, t + 5 * dt_6)
    f(k, tmp, p, t + 5 * dt_6)
    @.. broadcast=false thread=thread tmp=tmp+dt_6*k
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k, tmp, p, t + dt)
    @.. broadcast=false thread=thread u=(uprev+9*(u+dt_6*k₄)+
                                         15*(tmp+dt_6*k))/25
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 10)
end

function initialize!(integrator, cache::KYK2014DGSSPRK_3S2_ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return nothing
end

@muladd function perform_step!(integrator, cache::KYK2014DGSSPRK_3S2_ConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2) = cache
    u_1 = α_10 * uprev + dt * β_10 * integrator.fsalfirst
    u_2 = (α_20 * uprev +
           α_21 * u_1 + dt * β_21 * f(u_1, p, t + c_1 * dt))
    integrator.u = (α_30 * uprev + dt * β_30 * integrator.fsalfirst +
                    α_32 * u_2 + dt * β_32 * f(u_2, p, t + c_2 * dt))
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = f(integrator.u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.fsallast = integrator.k[2]
    return nothing
end

function initialize!(integrator, cache::KYK2014DGSSPRK_3S2_Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function perform_step!(integrator, cache::KYK2014DGSSPRK_3S2_Cache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, u_1, u_2, kk_1, kk_2, stage_limiter!, step_limiter!, thread) = cache
    (; α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2) = cache.tab

    @.. broadcast=false thread=thread u_1=α_10*uprev+dt*β_10*integrator.fsalfirst
    stage_limiter!(u_1, integrator, p, t + c_1 * dt)
    f(kk_1, u_1, p, t + c_1 * dt)
    @.. broadcast=false thread=thread u_2=(α_20*uprev+
                                           α_21*u_1+dt*β_21*kk_1)
    stage_limiter!(u_2, integrator, p, t + c_2 * dt)
    f(kk_2, u_2, p, t + c_2 * dt)
    @.. broadcast=false thread=thread u=(α_30*uprev+
                                         dt*β_30*integrator.fsalfirst+
                                         α_32*u_2+dt*β_32*kk_2)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(integrator.k[2], u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    return nothing
end
