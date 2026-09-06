function initialize!(integrator, cache::PDIRK44ConstantCache) end

@muladd function perform_step!(integrator, cache::PDIRK44ConstantCache, repeat_step = false)
    (; dt, uprev, u) = integrator
    alg = unwrap_alg(integrator, true)
    (; nlsolver, tab) = cache
    (; γs, cs, α1, α2, b1, b2, b3, b4) = tab

    k2 = Array{typeof(u)}(undef, 2)
    k1 = Array{typeof(u)}(undef, 2)
    let nlsolver = nlsolver, u = u, uprev = uprev, integrator = integrator,
            cache = cache, dt = dt, repeat_step = repeat_step,
            k1 = k1

        @threaded alg.threading for i in 1:2
            nlsolver[i].z = zero(u)
            nlsolver[i].tmp = uprev
            nlsolver[i].γ = γs[i]
            nlsolver[i].c = cs[i]
            markfirststage!(nlsolver[i])
            # The next stage may mutate the nonlinear solver's returned storage.
            k1[i] = copy(nlsolve!(nlsolver[i], integrator, cache, repeat_step))
        end
    end
    failed = any(nlsolvefail, nlsolver)
    merge_stats_deltas!(integrator.stats, nlsolver)
    integrator.force_stepfail = failed
    failed && return
    let nlsolver = nlsolver, u = u, uprev = uprev, integrator = integrator,
            cache = cache, dt = dt, repeat_step = repeat_step,
            k1 = k1, k2 = k2

        @threaded alg.threading for i in 1:2
            nlsolver[i].c = cs[2 + i]
            nlsolver[i].z = zero(u)
            nlsolver[i].tmp = uprev + α1[i] * k1[1] + α2[i] * k1[2]
            k2[i] = copy(nlsolve!(nlsolver[i], integrator, cache, repeat_step))
        end
    end
    failed = any(nlsolvefail, nlsolver)
    merge_stats_deltas!(integrator.stats, nlsolver)
    integrator.force_stepfail = failed
    failed && return
    integrator.u = uprev + b1 * k1[1] + b2 * k2[1] + b3 * k1[2] + b4 * k2[2]
end

function initialize!(integrator, cache::PDIRK44Cache) end

@muladd function perform_step!(integrator, cache::PDIRK44Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    alg = unwrap_alg(integrator, true)
    (; nlsolver, k1, k2, tab) = cache
    (; γs, cs, α1, α2, b1, b2, b3, b4) = tab
    let nlsolver = nlsolver, u = u, uprev = uprev, integrator = integrator,
            cache = cache, dt = dt, repeat_step = repeat_step,
            k1 = k1

        @threaded alg.threading for i in 1:2
            nlsolver[i].z .= zero(eltype(u))
            nlsolver[i].tmp .= uprev
            nlsolver[i].γ = γs[i]
            nlsolver[i].c = cs[i]
            markfirststage!(nlsolver[i])
            k1[i] .= nlsolve!(nlsolver[i], integrator, cache, repeat_step)
        end
    end
    failed = any(nlsolvefail, nlsolver)
    merge_stats_deltas!(integrator.stats, nlsolver)
    integrator.force_stepfail = failed
    failed && return
    let nlsolver = nlsolver, u = u, uprev = uprev, integrator = integrator,
            cache = cache, dt = dt, repeat_step = repeat_step,
            k1 = k1, k2 = k2

        @threaded alg.threading for i in 1:2
            nlsolver[i].c = cs[2 + i]
            nlsolver[i].z .= zero(eltype(u))
            @.. broadcast = false nlsolver[i].tmp = uprev + α1[i] * k1[1] + α2[i] * k1[2]
            k2[i] .= nlsolve!(nlsolver[i], integrator, cache, repeat_step)
        end
    end
    failed = any(nlsolvefail, nlsolver)
    merge_stats_deltas!(integrator.stats, nlsolver)
    integrator.force_stepfail = failed
    failed && return
    @.. broadcast = false u = uprev + b1 * k1[1] + b2 * k2[1] + b3 * k1[2] + b4 * k2[2]
end
