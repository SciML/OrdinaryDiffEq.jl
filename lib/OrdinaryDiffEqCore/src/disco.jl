function set_discontinuity(u, uprev, integrator)
    breakpointθ = find_discontinuity(u, uprev, integrator)
    dt = integrator.dt
    if 1.0e-10 < breakpointθ < 1.0
        return breakpointθ * dt
    end
    return -one(dt)
end

get_disco_probs(cache::AbstractControllerCache) = cache.controller.basic.disco_probs
get_disco_probs(cache::DummyControllerCache) = cache.disco_probs
get_disco_probs(cache::CompositeControllerCache) = get_disco_probs(first(cache.caches))

function find_discontinuity(u, uprev, integrator)
    cb = integrator.opts.callback
    dt = integrator.dt
    cb === nothing && return -one(dt)
    isempty(cb.continuous_callbacks) && return -one(dt)
    p = integrator.p
    t = integrator.t
    k = integrator.k
    breakpointθ = -one(dt)
    disco_probs = get_disco_probs(integrator.controller_cache)
    idx = 1
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
            _ode_addsteps!(disco_zero.k, disco_zero.tprev, disco_zero.uprev, disco_zero.u,
                disco_zero.dt, disco_zero.f, disco_zero.p, disco_zero.cache, false, true, false)
            for j in 1:len_cb
                if (disco_zero.out_low[j] * disco_zero.out_high[j] < zero(disco_zero.out_low[j]))
                    disco_zero.ind = j
                    sol = solve(disco_prob, tspan = (0.0, breakpointθ == -one(dt) ? 1.0 : breakpointθ); save_everystep = false)
                    tmp = sol[]
                    if (!isnan(tmp) && (breakpointθ < zero(breakpointθ) || tmp < breakpointθ))
                        breakpointθ = tmp
                        integrator.curr_discontinuity = idx
                        integrator.disco_checkpoint = integrator.t + tmp * dt
                    end
                end
            end
        else
            out_prev = i.condition(uprev, t, integrator)
            out_curr = i.condition(u, t + dt, integrator)
            if (out_prev * out_curr < zero(out_prev))
                _ode_addsteps!(disco_zero.k, disco_zero.tprev, disco_zero.uprev, disco_zero.u,
                            disco_zero.dt, disco_zero.f, disco_zero.p, disco_zero.cache, false, true, false)
                sol = solve(disco_prob, tspan = (0.0, breakpointθ == -one(dt) ? 1.0 : breakpointθ); save_everystep = false)
                tmp = sol[]
                if (!isnan(tmp) && (breakpointθ < zero(breakpointθ) || tmp < breakpointθ))
                    breakpointθ = tmp
                    integrator.curr_discontinuity = idx
                    integrator.disco_checkpoint = integrator.t + tmp * dt
                end
            end
        end
        idx += 1
    end
    return breakpointθ
end
