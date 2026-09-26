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

# The node rates only describe the step when it ends on the collocation polynomial.
sdc_stores_dense(integrator, alg) =
    integrator.opts.calck && alg.step_update === SDCStepUpdate.Quadrature

function initialize!(integrator, cache::SDCCache)
    sdc_stores_dense(integrator, unwrap_alg(integrator, true)) || return nothing
    integrator.kshortsize = 2
    resize!(integrator.k, 2)
    (; kdense) = cache
    integrator.f(kdense[1], integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    recursivecopy!(kdense[2], kdense[1])
    integrator.k[1] = kdense[1]
    integrator.k[2] = kdense[2]
    return nothing
end

function initialize!(integrator, cache::SDCConstantCache)
    sdc_stores_dense(integrator, unwrap_alg(integrator, true)) || return nothing
    integrator.kshortsize = 2
    du = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k = typeof(integrator.k)(undef, 2)
    integrator.k[1] = du
    integrator.k[2] = du
    return nothing
end

"""
    sdc_store_dense!(integrator, cache, zk, zEk, alg, M)

Hand the node rates of the step just taken to the interpolant, as `k[1]`, `k[2]` for the
derivatives at the step ends and `k[3:M + 2]` for the rates themselves.
"""
function sdc_store_dense!(integrator, cache::SDCCache, zk, zEk, alg, M)
    sdc_stores_dense(integrator, alg) || return nothing
    (; kdense) = cache
    k = integrator.k
    if length(k) != M + 2
        resize!(k, M + 2)
        for i in 1:(M + 2)
            k[i] = kdense[i]
        end
    end
    invdt = inv(integrator.dt)
    for m in 1:M
        km = kdense[m + 2]
        @.. broadcast = false km = invdt * zk[m]
        isempty(zEk) || @.. broadcast = false km = km + invdt * zEk[m]
    end
    (; dense) = cache.tab
    for (i, Θ) in ((1, false), (2, true))
        ki = kdense[i]
        fill!(ki, false)
        for m in 1:M
            w = sdc_basis_value(dense, m, Θ)
            km = kdense[m + 2]
            @.. broadcast = false ki = ki + w * km
        end
    end
    return nothing
end

function sdc_store_dense!(integrator, cache::SDCConstantCache, zk, zEk, alg, M)
    sdc_stores_dense(integrator, alg) || return nothing
    k = integrator.k
    length(k) == M + 2 || resize!(k, M + 2)
    invdt = inv(integrator.dt)
    for m in 1:M
        z = isempty(zEk) ? zk[m] : zk[m] + zEk[m]
        k[m + 2] = @.. broadcast = false invdt * z
    end
    (; dense) = cache.tab
    for i in 1:2
        Θ = i == 2
        w = sdc_basis_value(dense, 1, Θ)
        k3 = k[3]
        ki = @.. broadcast = false w * k3
        for m in 2:M
            w = sdc_basis_value(dense, m, Θ)
            km = k[m + 2]
            ki = @.. broadcast = false ki + w * km
        end
        k[i] = ki
    end
    return nothing
end

"""
    sdc_step_update!(u, uprev, weights, z, zE, ulast, step_update)

Form the step solution from the current node rates, in place. `zE` holds the
explicit rates of a split problem and is empty otherwise.
"""
function sdc_step_update!(u, uprev, weights, z, zE, ulast, step_update)
    if step_update === SDCStepUpdate.Quadrature
        @.. broadcast = false u = uprev
        for m in eachindex(weights)
            iszero(weights[m]) && continue
            @.. broadcast = false u = u + weights[m] * z[m]
            isempty(zE) || @.. broadcast = false u = u + weights[m] * zE[m]
        end
    else
        @.. broadcast = false u = ulast
    end
    return u
end

"""
    sdc_step_update(uprev, weights, z, zE, ulast, step_update)

Out-of-place counterpart of [`sdc_step_update!`](@ref).
"""
function sdc_step_update(uprev, weights, z, zE, ulast, step_update)
    step_update === SDCStepUpdate.Quadrature || return ulast
    u = uprev
    for m in eachindex(weights)
        iszero(weights[m]) && continue
        u = @.. broadcast = false u + weights[m] * z[m]
        isempty(zE) || (u = @.. broadcast = false u + weights[m] * zE[m])
    end
    return u
end

"""
    sdc_node!(m, integrator, cache, QΔ, zk, zk1, zEk, zEk1, repeat_step)

One node of one sweep, writing only into slot `m` of the per-node buffers.

Everything it reads is either shared and constant for the sweep (`zk`, `uprev`,
the coefficients) or private to node `m`. The one exception is the strictly
lower part of `QΔ`, which couples node `m` to nodes before it — that part is
empty for a diagonal `QΔ`, which is what makes the node loop safe to thread.
"""
@muladd function sdc_node!(m, integrator, cache, QΔ, zk, zk1, zEk, zEk1, repeat_step)
    (; t, dt, uprev, f, p) = integrator
    (; tmp, ubuf, k, k2, nlsolvers, tab, solver_index, split) = cache
    (; nodes, Q, QE) = tab
    M = length(nodes)
    tm = t + nodes[m] * dt

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
    if split
        for j in 1:M
            coeff = Q[m, j] - QE[m, j]
            iszero(coeff) && continue
            @.. broadcast = false tmpm = tmpm + coeff * zEk[j]
        end
        for j in 1:(m - 1)
            coeff = QE[m, j]
            iszero(coeff) && continue
            @.. broadcast = false tmpm = tmpm + coeff * zEk1[j]
        end
    end

    index = solver_index[m]
    if iszero(index)
        # QΔ[m,m] = 0, so the node is explicit and u_m is the right-hand side.
        @.. broadcast = false ubuf[m] = tmpm
        split ? f.f1(k[m], ubuf[m], p, tm) : f(k[m], ubuf[m], p, tm)
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
    if split
        f.f2(k2[m], ubuf[m], p, tm)
        cache.nf2[m] += 1
        @.. broadcast = false zEk1[m] = dt * k2[m]
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::SDCCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; ubuf, ulow, atmp, k, k2, tab, failed, split) = cache
    (; nodes, weights) = tab
    alg = unwrap_alg(integrator, true)
    M = length(nodes)
    threading = alg.threading

    # COPY initialisation: u⁰_m = u_n at every node, which is what the standard
    # order predictions for SDC assume.
    for m in 1:M
        tm = t + nodes[m] * dt
        if split
            f.f1(k[m], uprev, p, tm)
            f.f2(k2[m], uprev, p, tm)
            @.. broadcast = false cache.zE[m] = dt * k2[m]
        else
            f(k[m], uprev, p, tm)
        end
        @.. broadcast = false cache.z[m] = dt * k[m]
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, M)
    split && (integrator.stats.nf2 += M)
    @.. broadcast = false ubuf[M] = uprev

    zk, zk1 = cache.z, cache.znew
    zEk, zEk1 = cache.zE, cache.zE_new
    adaptive = integrator.opts.adaptive
    # The step update after sweep k-1 is the embedded solution, so it is formed
    # every sweep and kept one behind.
    adaptive && sdc_step_update!(u, uprev, weights, zk, zEk, ubuf[M], alg.step_update)
    for sweep in 1:(alg.num_sweeps)
        QΔ = sdc_qdelta_for(tab, sweep)
        fill!(failed, false)
        fill!(cache.nf, 0)
        fill!(cache.nf2, 0)
        # `let` so the closure the threading macro builds captures the current
        # sweep's arrays by value rather than boxing the rebound names.
        let cache = cache, integrator = integrator, QΔ = QΔ, zk = zk, zk1 = zk1,
                zEk = zEk, zEk1 = zEk1, repeat_step = repeat_step

            @threaded threading for m in 1:M
                sdc_node!(m, integrator, cache, QΔ, zk, zk1, zEk, zEk1, repeat_step)
            end
        end
        # `nlsolve!` writes `integrator.force_stepfail` from every node, so a
        # later success would erase an earlier failure. The per-node flags are
        # the reliable record.
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, sum(cache.nf))
        split && (integrator.stats.nf2 += sum(cache.nf2))
        if any(failed)
            integrator.force_stepfail = true
            return nothing
        end
        zk, zk1 = zk1, zk
        zEk, zEk1 = zEk1, zEk
        adaptive && @.. broadcast = false ulow = u
        adaptive && sdc_step_update!(u, uprev, weights, zk, zEk, ubuf[M], alg.step_update)
    end

    adaptive || sdc_step_update!(u, uprev, weights, zk, zEk, ubuf[M], alg.step_update)
    sdc_store_dense!(integrator, cache, zk, zEk, alg, M)

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
    (; nlsolvers, tab, solver_index, split) = cache
    (; nodes, weights, Q, QE) = tab
    alg = unwrap_alg(integrator, true)
    M = length(nodes)
    stats = integrator.stats

    fimpl = split ? f.f1 : f
    zk = [dt * fimpl(uprev, p, t + nodes[m] * dt) for m in 1:M]
    zk1 = copy(zk)
    OrdinaryDiffEqCore.increment_nf!(stats, M)
    zEk = split ? [dt * f.f2(uprev, p, t + nodes[m] * dt) for m in 1:M] : empty(zk)
    zEk1 = copy(zEk)
    split && (stats.nf2 += M)
    ulast = uprev

    adaptive = integrator.opts.adaptive
    u = adaptive ? sdc_step_update(uprev, weights, zk, zEk, ulast, alg.step_update) : uprev
    ulow = u
    for sweep in 1:(alg.num_sweeps)
        QΔ = sdc_qdelta_for(tab, sweep)
        for m in 1:M
            tm = t + nodes[m] * dt
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
            if split
                for j in 1:M
                    coeff = Q[m, j] - QE[m, j]
                    iszero(coeff) && continue
                    tmp = @.. broadcast = false tmp + coeff * zEk[j]
                end
                for j in 1:(m - 1)
                    coeff = QE[m, j]
                    iszero(coeff) && continue
                    tmp = @.. broadcast = false tmp + coeff * zEk1[j]
                end
            end
            index = solver_index[m]
            if iszero(index)
                ulast = tmp
                zk1[m] = dt * fimpl(ulast, p, tm)
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
            if split
                zEk1[m] = dt * f.f2(ulast, p, tm)
                stats.nf2 += 1
            end
        end
        zk, zk1 = zk1, zk
        zEk, zEk1 = zEk1, zEk
        if adaptive
            ulow = u
            u = sdc_step_update(uprev, weights, zk, zEk, ulast, alg.step_update)
        end
    end

    adaptive || (u = sdc_step_update(uprev, weights, zk, zEk, ulast, alg.step_update))
    integrator.u = u
    sdc_store_dense!(integrator, cache, zk, zEk, alg, M)

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
