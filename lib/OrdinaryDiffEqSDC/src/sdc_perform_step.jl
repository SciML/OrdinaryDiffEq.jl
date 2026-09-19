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

function initialize!(integrator, cache::SDCCache)
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
    integrator.kshortsize = 2
    du = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k = typeof(integrator.k)(undef, 2)
    integrator.k[1] = du
    integrator.k[2] = du
    return nothing
end

"""
    sdc_store_dense!(integrator, cache, zk, alg, M)

Hand the node rates of the step just taken to the interpolant, as `k[1]`, `k[2]` for the
derivatives at the step ends and `k[3:M + 2]` for the rates themselves.

Only the quadrature update ends the step on the collocation polynomial, so `LastNode` keeps
the generic interpolation, and nothing is stored when no interpolation was asked for.
"""
function sdc_store_dense!(integrator, cache::SDCCache, zk, alg, M)
    (integrator.opts.calck && alg.step_update === SDCStepUpdate.Quadrature) ||
        return nothing
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
        @.. broadcast = false kdense[m + 2] = invdt * zk[m]
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

function sdc_store_dense!(integrator, cache::SDCConstantCache, zk, alg, M)
    (integrator.opts.calck && alg.step_update === SDCStepUpdate.Quadrature) ||
        return nothing
    k = integrator.k
    length(k) == M + 2 || resize!(k, M + 2)
    invdt = inv(integrator.dt)
    for m in 1:M
        k[m + 2] = @.. broadcast = false invdt * zk[m]
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

@muladd function perform_step!(integrator, cache::SDCCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, ubuf, ulow, atmp, k, nlsolvers, tab, solver_index) = cache
    (; nodes, weights, Q) = tab
    alg = unwrap_alg(integrator, true)
    M = length(nodes)
    stats = integrator.stats

    # COPY initialisation: u⁰_m = u_n at every node, which is what the standard
    # order predictions for SDC assume.
    for m in 1:M
        f(k, uprev, p, t + nodes[m] * dt)
        @.. broadcast = false cache.z[m] = dt * k
    end
    OrdinaryDiffEqCore.increment_nf!(stats, M)
    @.. broadcast = false ubuf = uprev

    zk, zk1 = cache.z, cache.znew
    adaptive = integrator.opts.adaptive
    # The step update after sweep k-1 is the embedded solution, so it is formed
    # every sweep and kept one behind.
    adaptive && sdc_step_update!(u, uprev, weights, zk, ubuf, alg.step_update)
    for sweep in 1:(alg.num_sweeps)
        QΔ = sdc_qdelta_for(tab, sweep)
        for m in 1:M
            @.. broadcast = false tmp = uprev
            for j in 1:M
                coeff = Q[m, j] - QΔ[m, j]
                iszero(coeff) && continue
                @.. broadcast = false tmp = tmp + coeff * zk[j]
            end
            for j in 1:(m - 1)
                coeff = QΔ[m, j]
                iszero(coeff) && continue
                @.. broadcast = false tmp = tmp + coeff * zk1[j]
            end
            index = solver_index[m]
            if iszero(index)
                # QΔ[m,m] = 0, so the node is explicit and u_m is the right-hand side.
                @.. broadcast = false ubuf = tmp
                f(k, ubuf, p, t + nodes[m] * dt)
                OrdinaryDiffEqCore.increment_nf!(stats, 1)
                @.. broadcast = false zk1[m] = dt * k
            else
                nls = nlsolvers[index]
                @.. broadcast = false nls.tmp = tmp
                @.. broadcast = false nls.z = zk[m]
                nls.γ = QΔ[m, m]
                nls.c = nodes[m]
                markfirststage!(nls)
                znode = nlsolve!(nls, integrator, cache, repeat_step)
                nlsolvefail(nls) && return
                @.. broadcast = false zk1[m] = znode
                @.. broadcast = false ubuf = tmp + QΔ[m, m] * znode
            end
        end
        zk, zk1 = zk1, zk
        adaptive && @.. broadcast = false ulow = u
        adaptive && sdc_step_update!(u, uprev, weights, zk, ubuf, alg.step_update)
    end

    adaptive || sdc_step_update!(u, uprev, weights, zk, ubuf, alg.step_update)
    sdc_store_dense!(integrator, cache, zk, alg, M)

    if adaptive
        @.. broadcast = false tmp = u - ulow
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
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
    sdc_store_dense!(integrator, cache, zk, alg, M)

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
