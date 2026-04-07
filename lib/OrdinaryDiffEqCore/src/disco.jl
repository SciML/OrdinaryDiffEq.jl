function set_discontinuity(u, uprev, integrator, cache) 
    breakpointθ = find_discontinuity(u, uprev, integrator, cache) 
    dt = integrator.dt
    t = integrator.t
    if !isnan(breakpointθ) && 1e-6 < breakpointθ < 1.0
        #println("Discontinuity detected at t = ", t + breakpointθ * dt)
        return breakpointθ * dt
    end
    return -1
end

function find_discontinuity(u, uprev, integrator, cache)
    cb = integrator.opts.callback
    cb === nothing && return -1
    isempty(cb.continuous_callbacks) && return -1
    p = integrator.p
    t = integrator.t
    dt = integrator.dt
    save_idxs = integrator.opts.save_idxs
    k = integrator.k
    cache = integrator.cache
    differential_vars = integrator.differential_vars
    θlo = zero(dt)
    θhi = one(dt)
    bracket = [θlo, θhi]
    breakpointθ = -one(dt)
    idx = 1
    for i in cb.continuous_callbacks
        if (!(i.is_discontinuity)) 
            continue 
        end
        disco_prob = integrator.disco_probs[idx]
        disco_zero = disco_prob.f.f
        disco_zero.dt = dt
        disco_zero.uprev = uprev
        disco_zero.u = u
        disco_zero.k = k
        disco_zero.cache = cache
        disco_zero.differential_vars = differential_vars
        disco_zero.idxs = save_idxs
        if (i isa VectorContinuousCallback)
            len_cb = i.len
            out_prev = similar(u)
            out_curr = similar(u)
            i.condition(out_prev, uprev, t, integrator)
            i.condition(out_curr, u, t + dt, integrator)
            for j in 1:len_cb
                if (out_prev[j] * out_curr[j] < zero(out_prev[j]))
                    disco_zero.ind = j
                    sol = solve(disco_prob; bracket = bracket)
                    tmp = sol[]
                    if (!isnan(tmp) && (breakpointθ == -1 || tmp < breakpointθ)) 
                        breakpointθ = tmp 
                    end 
                end
            end
        else
            out_prev = i.condition(uprev, t, integrator)
            out_curr = i.condition(u, t + dt, integrator)
            if (out_prev * out_curr < zero(out_prev))
                sol = solve(disco_prob; bracket = bracket)
                tmp = sol[]
                if (!isnan(tmp) && (breakpointθ == -1 || tmp < breakpointθ)) 
                    breakpointθ = tmp 
                end 
            end
        end
        idx += 1
    end
    breakpointθ
end
