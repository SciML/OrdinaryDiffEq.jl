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
    breakpointθ = -one(dt)
    idx = 1
    for i in cb.continuous_callbacks
        if (!(i.is_discontinuity)) 
            continue 
        end
        if (i isa VectorContinuousCallback)
            out_prev = similar(u)
            out_curr = similar(u)  
            i.condition(out_prev, uprev, t, integrator)
            i.condition(out_curr, u, t + dt, integrator)
            for (ind, (f0, f1)) in enumerate(zip(out_prev, out_curr))
                if (f0 * f1 < zero(f0))
                    u₁ = similar(u)
                    out = similar(u)
                    function zero_func(θ, p)
                        ode_interpolant!(u₁, θ, integrator, integrator.opts.save_idxs, Val{0})
                        i.condition(out, u₁, t + θ * integrator.dt, integrator)
                        out[ind]
                    end
                    prob = IntervalNonlinearProblem(zero_func, [zero(dt), one(dt)], p)
                    sol = solve(prob; bracket=[zero(dt), one(dt)], abstol = 0, reltol = 0)                
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
                disco_prob = integrator.disco_probs[idx]
                #disco_prob = integrator.disco_prob
                disco_prob.f.f.dt = integrator.dt
                disco_prob.f.f.uprev = uprev
                disco_prob.f.f.u = u
                disco_prob.f.f.k = integrator.k
                disco_prob.f.f.cache = integrator.cache
                disco_prob.f.f.differential_vars = integrator.differential_vars
                disco_prob.f.f.idxs = integrator.opts.save_idxs
                #disco_prob.f.f.callback = i                
                sol = solve(disco_prob; bracket=[zero(dt), one(dt)], abstol = 0, reltol = 0)                
                tmp = sol[]
                if (!isnan(tmp) && (breakpointθ == -1 || tmp < breakpointθ)) 
                    breakpointθ = tmp 
                end 
            end
            idx += 1
        end
    end
    breakpointθ
end
