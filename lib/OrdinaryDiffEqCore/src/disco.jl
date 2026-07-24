"""
    set_discontinuity(integrator) -> dt

Return the sub-step `Δt` at which a continuous callback's root lies within the
current step (a discontinuity to stop at), or a negative value if there is none in
`(0, 1)·dt`.
"""
function set_discontinuity(integrator)
    breakpointθ = find_discontinuity(integrator)
    dt = integrator.dt
    if breakpointθ < one(breakpointθ)
        # below formula is equivalent to (0.5 + 0.4 * sin(π * (breakpointθ - 0.5))) * dt
        return (0.5 + 1.2372 * (breakpointθ - 0.5) - 1.7487 * (breakpointθ - 0.5)^3) * dt
    end
    return -one(dt)
end

get_disco_probs(cache::AbstractControllerCache) = cache.controller.basic.disco_probs
get_disco_probs(cache::DummyControllerCache) = cache.disco_probs
get_disco_probs(cache::CompositeControllerCache) = get_disco_probs(first(cache.caches))

function find_discontinuity(integrator)
    cb = integrator.opts.callback
    dt = integrator.dt
    u = integrator.u
    uprev = integrator.uprev
    cb === nothing && return one(dt)
    isempty(cb.continuous_callbacks) && return one(dt)
    p = integrator.p
    t = integrator.t
    k = integrator.k
    breakpointθ = one(dt)
    disco_probs = get_disco_probs(integrator.controller_cache)
    idx = 1
    addsteps_called = false
    for i in cb.continuous_callbacks
        if (!(i.maybe_discontinuity))
            continue
        end
        disco_prob = disco_probs[idx]
        disco_zero = disco_prob.f.f.obj.x
        disco_zero.dt = dt
        disco_zero.uprev = uprev
        disco_zero.u = u
        disco_zero.k = k
        disco_zero.tprev = t
        disco_zero.p = p
        if (i isa VectorContinuousCallback)
            len_cb = i.len
            i.condition(disco_zero.out_low, uprev, t, integrator)
            i.condition(disco_zero.out_high, u, t + dt, integrator)
            for j in 1:len_cb
                if (disco_zero.out_low[j] * disco_zero.out_high[j] < zero(disco_zero.out_low[j]))
                    if (!addsteps_called)
                        addsteps_called = true
                        _ode_addsteps!(disco_zero.k, disco_zero.tprev, disco_zero.uprev, disco_zero.u,
                                    disco_zero.dt, disco_zero.f, disco_zero.p, disco_zero.cache, false, true, false)
                    end
                    disco_zero.ind = j
                    disco_prob.tspan[2] = breakpointθ
                    sol = solve(disco_prob)
                    tmp = sol[]
                    if (!isnan(tmp) && tmp < breakpointθ)
                        breakpointθ = tmp
                        integrator.is_disco_step = true
                        integrator.disco_checkpoint = integrator.t + dt #our prev rejected step, we shouldn't step too far past this
                    end
                end
            end
        else
            disco_zero.out_low[1] = i.condition(uprev, t, integrator)
            disco_zero.out_high[1] = i.condition(u, t + dt, integrator)
            if (disco_zero.out_low[1] * disco_zero.out_high[1] < zero(disco_zero.out_low[1]))
                if (!addsteps_called)
                    addsteps_called = true
                    _ode_addsteps!(disco_zero.k, disco_zero.tprev, disco_zero.uprev, disco_zero.u,
                                disco_zero.dt, disco_zero.f, disco_zero.p, disco_zero.cache, false, true, false)
                end
                disco_prob.tspan[2] = breakpointθ
                sol = solve(disco_prob)
                tmp = sol[]
                if (!isnan(tmp) && tmp < breakpointθ)
                    breakpointθ = tmp
                    integrator.is_disco_step = true
                    integrator.disco_checkpoint = integrator.t + dt #our prev rejected step, we shouldn't step past this
                end
            end
        end
        idx += 1
    end
    return breakpointθ
end
