# The sweep implemented below is
#
#     u_m^{k+1} = rhs_m + Δt QΔ[m,m] f(u_m^{k+1}, t_n + τ_m Δt),
#     rhs_m = u_n + Δt Σ_{j=1}^{M} (Q - QΔ)[m,j] f_j^k
#                 + Δt Σ_{j<m}     QΔ[m,j]       f_j^{k+1}.
#
# `z_m = Δt f_m` is stored rather than `f_m` or `u_m` because that is the
# variable `nlsolve!` already solves for.
#
# The `Σ_{j<m}` term is the only node-to-node coupling, and it is empty for a
# diagonal `QΔ`. Keeping the loop body dependent on nothing but `z` (sweep `k`),
# `znew[1:m-1]` and node-local scratch is what will let the `m` loop be threaded
# for parallel-across-the-nodes SDC without restructuring it.

function initialize!(integrator, cache::Union{SDCCache, SDCConstantCache}) end

"""
    sdc_step_update!(u, uprev, weights, z, ulast, step_update)

Form the step solution from the current node rates, in place.
"""
function sdc_step_update!(u, uprev, weights, z, ulast, step_update)
    if step_update === SDCStepUpdate.Quadrature
        @.. broadcast = false u = uprev
        for m in eachindex(weights)
            iszero(weights[m]) && continue
            @.. broadcast = false u = u + weights[m] * z[m]
        end
    else
        @.. broadcast = false u = ulast
    end
    return u
end

"""
    sdc_step_update(uprev, weights, z, ulast, step_update)

Out-of-place counterpart of [`sdc_step_update!`](@ref).
"""
function sdc_step_update(uprev, weights, z, ulast, step_update)
    step_update === SDCStepUpdate.Quadrature || return ulast
    u = uprev
    for m in eachindex(weights)
        iszero(weights[m]) && continue
        u = @.. broadcast = false u + weights[m] * z[m]
    end
    return u
end

"""
    sdc_node!(m, integrator, cache, QΔ, zk, zk1, repeat_step)

One node of one sweep, writing only into slot `m` of the per-node buffers.

Everything it reads is either shared and constant for the sweep (`zk`, `uprev`,
the coefficients) or private to node `m`. The one exception is the strictly
lower part of `QΔ`, which couples node `m` to nodes before it — that part is
empty for a diagonal `QΔ`, which is what makes the node loop safe to thread.
"""
@muladd function sdc_node!(m, integrator, cache, QΔ, zk, zk1, repeat_step)
    (; t, dt, uprev, f, p) = integrator
    (; tmp, ubuf, k, nlsolvers, tab, solver_index) = cache
    (; nodes, Q) = tab
    M = length(nodes)

    tmpm = tmp[m]
    @.. broadcast = false tmpm = uprev
    for j in 1:M
        coeff = Q[m, j] - QΔ[m, j]
        iszero(coeff) && continue
        @.. broadcast = false tmpm = tmpm + coeff * zk[j]
    end
    for j in 1:(m - 1)
        coeff = QΔ[m, j]
        iszero(coeff) && continue
        @.. broadcast = false tmpm = tmpm + coeff * zk1[j]
    end

    index = solver_index[m]
    if iszero(index)
        # QΔ[m,m] = 0, so the node is explicit and u_m is the right-hand side.
        @.. broadcast = false ubuf[m] = tmpm
        f(k[m], ubuf[m], p, t + nodes[m] * dt)
        cache.nf[m] += 1
        @.. broadcast = false zk1[m] = dt * k[m]
    else
        nls = nlsolvers[index]
        @.. broadcast = false nls.tmp = tmpm
        @.. broadcast = false nls.z = zk[m]
        nls.γ = QΔ[m, m]
        nls.c = nodes[m]
        markfirststage!(nls)
        znode = nlsolve!(nls, integrator, cache, repeat_step)
        cache.failed[m] = nlsolvefail(nls)
        @.. broadcast = false zk1[m] = znode
        @.. broadcast = false ubuf[m] = tmpm + QΔ[m, m] * znode
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::SDCCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; ubuf, ulow, atmp, k, tab, failed) = cache
    (; nodes, weights) = tab
    alg = unwrap_alg(integrator, true)
    M = length(nodes)
    threading = alg.threading

    for m in 1:M
        f(k[m], uprev, p, t + nodes[m] * dt)
        @.. broadcast = false cache.z[m] = dt * k[m]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, M)
    @.. broadcast = false ubuf[M] = uprev

    zk, zk1 = cache.z, cache.znew
    adaptive = integrator.opts.adaptive
    adaptive && sdc_step_update!(u, uprev, weights, zk, ubuf[M], alg.step_update)
    for sweep in 1:(alg.num_sweeps)
        QΔ = sdc_qdelta_for(tab, sweep)
        fill!(failed, false)
        fill!(cache.nf, 0)
        # `let` so the closure the threading macro builds captures the current
        # sweep's arrays by value rather than boxing the rebound names.
        let cache = cache, integrator = integrator, QΔ = QΔ, zk = zk, zk1 = zk1,
                repeat_step = repeat_step

            @threaded threading for m in 1:M
                sdc_node!(m, integrator, cache, QΔ, zk, zk1, repeat_step)
            end
        end
        # `nlsolve!` writes `integrator.force_stepfail` from every node, so a
        # later success would erase an earlier failure. The per-node flags are
        # the reliable record.
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, sum(cache.nf))
        if any(failed)
            integrator.force_stepfail = true
            return nothing
        end
        zk, zk1 = zk1, zk
        adaptive && @.. broadcast = false ulow = u
        adaptive && sdc_step_update!(u, uprev, weights, zk, ubuf[M], alg.step_update)
    end

    adaptive || sdc_step_update!(u, uprev, weights, zk, ubuf[M], alg.step_update)

    if adaptive
        tmp1 = cache.tmp[1]
        @.. broadcast = false tmp1 = u - ulow
        calculate_residuals!(
            atmp, tmp1, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::SDCConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    (; nlsolvers, tab, solver_index) = cache
    (; nodes, weights, Q) = tab
    alg = unwrap_alg(integrator, true)
    M = length(nodes)
    stats = integrator.stats

    zk = [dt * f(uprev, p, t + nodes[m] * dt) for m in 1:M]
    zk1 = copy(zk)
    OrdinaryDiffEqCore.increment_nf!(stats, M)
    ulast = uprev

    adaptive = integrator.opts.adaptive
    u = adaptive ? sdc_step_update(uprev, weights, zk, ulast, alg.step_update) : uprev
    ulow = u
    for sweep in 1:(alg.num_sweeps)
        QΔ = sdc_qdelta_for(tab, sweep)
        for m in 1:M
            tmp = uprev
            for j in 1:M
                coeff = Q[m, j] - QΔ[m, j]
                iszero(coeff) && continue
                tmp = @.. broadcast = false tmp + coeff * zk[j]
            end
            for j in 1:(m - 1)
                coeff = QΔ[m, j]
                iszero(coeff) && continue
                tmp = @.. broadcast = false tmp + coeff * zk1[j]
            end
            index = solver_index[m]
            if iszero(index)
                ulast = tmp
                zk1[m] = dt * f(ulast, p, t + nodes[m] * dt)
                OrdinaryDiffEqCore.increment_nf!(stats, 1)
            else
                nls = nlsolvers[index]
                nls.tmp = tmp
                nls.z = zk[m]
                nls.γ = QΔ[m, m]
                nls.c = nodes[m]
                markfirststage!(nls)
                znode = nlsolve!(nls, integrator, cache, repeat_step)
                nlsolvefail(nls) && return
                zk1[m] = znode
                ulast = @.. broadcast = false tmp + QΔ[m, m] * znode
            end
        end
        zk, zk1 = zk1, zk
        if adaptive
            ulow = u
            u = sdc_step_update(uprev, weights, zk, ulast, alg.step_update)
        end
    end

    adaptive || (u = sdc_step_update(uprev, weights, zk, ulast, alg.step_update))
    integrator.u = u

    if adaptive
        utilde = @.. broadcast = false u - ulow
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    return nothing
end
