function set_discontinuity(u, uprev, integrator, cache) #need to pick algs to test
    breakpointθ = find_discontinuity(u, uprev, integrator, cache) 
    dt = integrator.dt
    t = integrator.t
    if !isnan(breakpointθ) && 1e-6 < breakpointθ < 1.0
        #println("Discontinuity detected at t = ", t + breakpointθ * dt)
        integrator.dt = breakpointθ * dt
        integrator.disco_dt_set = true
    end
end

function find_discontinuity(u, uprev, integrator, cache)
    cb = integrator.opts.callback
    cb === nothing && return -1
    isempty(cb.continuous_callbacks) && return -1

    disco_exists = false;
    for i in cb.continuous_callbacks
        if (i.is_discontinuity) 
            disco_exists = true
            break
        end
    end
    !disco_exists && return -1
    p = integrator.p
    t = integrator.t
    dt = integrator.dt
    breakpointθ = -one(dt)
    prob = nothing
    for i in cb.continuous_callbacks
        if (!(i.is_discontinuity)) 
            continue 
        end
        out_prev = nothing
        out_curr = nothing
        is_inplace = DiffEqBase.isinplace(i.condition, 4)
        if is_inplace
            out_prev = similar(u)
            i.condition(out_prev, uprev, t, integrator)
            out_curr = similar(u)
            i.condition(out_curr, u, t + dt, integrator)
            is_inplace = true
        else
            out_prev = i.condition(uprev, t, integrator)
            out_curr = i.condition(u, t + dt, integrator)
            is_inplace = false
        end
        for (idx, (f0, f1)) in enumerate(zip(out_prev, out_curr))
            if (f0 * f1 < zero(f0))
                function zero_func(θ, p)
                    u₁ = similar(u)
                    _ode_interpolant!(u₁, θ, dt, uprev, u, integrator.k, cache,
                                    nothing, Val{0}, nothing)

                    if is_inplace
                        out = similar(u)
                        i.condition(out, u₁, t + θ * dt, integrator)
                    else
                        out = i.condition(u₁, t + θ * dt, integrator)
                    end
                    out[idx]
                end
                if prob === nothing
                    prob = IntervalNonlinearProblem(zero_func, [zero(dt), one(dt)], p)
                else
                    prob = remake(prob; f=zero_func)
                end
                sol = solve(prob; bracket=[zero(dt), one(dt)])                
                tmp = sol[]
                if (!isnan(tmp) && (breakpointθ == -1 || tmp < breakpointθ)) 
                    breakpointθ = tmp 
                end
            end
        end
    end
    breakpointθ
end

