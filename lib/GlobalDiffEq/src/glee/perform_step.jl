function initialize!(integrator, cache::GLEECache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.fsalfirst
    integrator.k[2] = cache.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function perform_step!(integrator, cache::GLEECache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; ks, ytmp, εloc, atmp, tab) = cache
    inner = _glee_inner(f)
    yprev = uprev.x[1]
    εprev = uprev.x[2]
    s = nstages(tab)

    for i in 1:s
        if i == 1 && tab.stage1_fsal
            copyto!(ks[1], integrator.fsalfirst.x[1])
            continue
        end
        @.. ytmp = tab.U[i, 1] * yprev + tab.U[i, 2] * εprev
        for j in 1:(i - 1)
            iszero(tab.A[i, j]) && continue
            @.. ytmp = ytmp + dt * tab.A[i, j] * ks[j]
        end
        inner(ks[i], ytmp, p, t + tab.c[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    ynew = u.x[1]
    εnew = u.x[2]
    fill!(εloc, zero(eltype(εloc)))
    for j in 1:s
        iszero(tab.B[2, j]) && continue
        @.. εloc = εloc + dt * tab.B[2, j] * ks[j]
    end
    copyto!(ynew, yprev)
    for j in 1:s
        iszero(tab.B[1, j]) && continue
        @.. ynew = ynew + dt * tab.B[1, j] * ks[j]
    end
    @.. εnew = εprev + εloc

    if tab.solution_stage == 0
        f(integrator.fsallast, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    else
        copyto!(integrator.fsallast.x[1], ks[tab.solution_stage])
    end
    εslope = integrator.fsallast.x[2]
    @.. εslope = εloc / dt

    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, εloc, yprev, ynew, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    return nothing
end

function initialize!(integrator, cache::GLEEConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

@muladd function perform_step!(integrator, cache::GLEEConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    tab = cache.tab
    inner = _glee_inner(f)
    yprev = uprev.x[1]
    εprev = uprev.x[2]
    s = nstages(tab)

    k1 = if tab.stage1_fsal
        integrator.fsalfirst.x[1]
    else
        Y1 = tab.U[1, 1] * yprev + tab.U[1, 2] * εprev
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        inner(Y1, p, t + tab.c[1] * dt)
    end
    ks = Vector{typeof(k1)}(undef, s)
    ks[1] = k1
    for i in 2:s
        Y = tab.U[i, 1] * yprev + tab.U[i, 2] * εprev
        for j in 1:(i - 1)
            iszero(tab.A[i, j]) && continue
            Y = Y + dt * tab.A[i, j] * ks[j]
        end
        ks[i] = inner(Y, p, t + tab.c[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    ynew = yprev
    εloc = zero(yprev)
    for j in 1:s
        ynew = ynew + dt * tab.B[1, j] * ks[j]
        εloc = εloc + dt * tab.B[2, j] * ks[j]
    end
    εnew = εprev + εloc
    integrator.u = RecursiveArrayTools.ArrayPartition(ynew, εnew)

    fy = if tab.solution_stage == 0
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        inner(ynew, p, t + dt)
    else
        ks[tab.solution_stage]
    end
    integrator.fsallast = RecursiveArrayTools.ArrayPartition(fy, εloc / dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        atmp = calculate_residuals(
            εloc, yprev, ynew, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    return nothing
end
