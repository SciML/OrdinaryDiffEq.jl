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

@muladd function perform_step!(integrator, cache::SDCCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, ubuf, k, nlsolvers, tab, solver_index) = cache
    (; nodes, weights, Q, QΔ) = tab
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
    for _ in 1:(alg.num_sweeps)
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
    end

    if alg.step_update === SDCStepUpdate.Quadrature
        @.. broadcast = false u = uprev
        for m in 1:M
            iszero(weights[m]) && continue
            @.. broadcast = false u = u + weights[m] * zk[m]
        end
    else
        @.. broadcast = false u = ubuf
    end
    return nothing
end

@muladd function perform_step!(integrator, cache::SDCConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    (; nlsolvers, tab, solver_index) = cache
    (; nodes, weights, Q, QΔ) = tab
    alg = unwrap_alg(integrator, true)
    M = length(nodes)
    stats = integrator.stats

    zk = [dt * f(uprev, p, t + nodes[m] * dt) for m in 1:M]
    zk1 = copy(zk)
    OrdinaryDiffEqCore.increment_nf!(stats, M)
    ulast = uprev

    for _ in 1:(alg.num_sweeps)
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
    end

    integrator.u = if alg.step_update === SDCStepUpdate.Quadrature
        unew = uprev
        for m in 1:M
            iszero(weights[m]) && continue
            unew = @.. broadcast = false unew + weights[m] * zk[m]
        end
        unew
    else
        ulast
    end
    return nothing
end
