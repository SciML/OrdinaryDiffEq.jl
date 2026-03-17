"""
    track_propagated_discontinuities!(integrator::DDEIntegrator)

Try to find a propagated discontinuity in the time interval `[integrator.t, integrator.t +
integrator.dt]` and add it to the set of discontinuities and grid points of the
`integrator`.
"""
function track_propagated_discontinuities!(integrator::DDEIntegrator)
    # calculate interpolation points
    interp_points = integrator.discontinuity_interp_points
    Θs = range(zero(integrator.t); stop = oneunit(integrator.t), length = interp_points)

    @SciMLMessage(
        lazy"Tracking discontinuities for state-dependent delays in interval [$(integrator.t), $(integrator.t + integrator.dt)]",
        integrator.opts.verbose, :state_dependent_delay
    )

    # for dependent lags and previous discontinuities
    for lag in integrator.sol.prob.dependent_lags,
            discontinuity in integrator.tracked_discontinuities
        # obtain time of previous discontinuity
        T = discontinuity.t

        # estimate subinterval of current integration step that contains a propagated
        # discontinuity induced by the lag and the previous discontinuity
        interval = discontinuity_interval(integrator, lag, T, Θs)

        # if a discontinuity exists in the current integration step
        if interval !== nothing
            # estimate time point of discontinuity
            t = discontinuity_time(integrator, lag, T, interval)

            # add new discontinuity of correct order at the estimated time point
            if integrator.sol.prob.neutral
                d = Discontinuity(t, discontinuity.order)
                @SciMLMessage(
                    lazy"Propagated discontinuity found at t = $t with order $(discontinuity.order) (neutral DDE)",
                    integrator.opts.verbose, :discontinuity_tracking
                )
            else
                d = Discontinuity(t, discontinuity.order + 1)
                @SciMLMessage(
                    lazy"Propagated discontinuity found at t = $t with order $(discontinuity.order + 1)",
                    integrator.opts.verbose, :discontinuity_tracking
                )
            end
            push!(integrator.opts.d_discontinuities, d)
            push!(integrator.opts.tstops, t)

            # analogously to RADAR5 we do not strive for finding the first discontinuity
            break
        end
    end

    return nothing
end

"""
    discontinuity_function(integrator::DDEIntegrator, lag, T, t)

Evaluate function ``f(x) = T + lag(u(x), p, x) - x`` at time point `t`, where `T` is time
point of a previous discontinuity and `lag` is a dependent delay.
"""
function discontinuity_function(integrator::DDEIntegrator, lag, T, t)
    tmp = get_tmp_cache(integrator)
    cache = tmp === nothing ? nothing : first(tmp)

    # estimate state at the given time
    if cache === nothing
        ut = integrator(t, Val{0})
    else
        integrator(cache, t, Val{0})
        ut = cache
    end

    lag_value = lag(ut, integrator.p, t)
    @SciMLMessage(
        lazy"Evaluating state-dependent delay at t = $t: lag = $lag_value",
        integrator.opts.verbose, :delay_evaluation
    )
    return T + lag_value - t
end

"""
    discontinuity_interval(integrator::DDEIntegrator, lag, T, Θs)

Return an estimated subinterval of the current integration step of the `integrator` that
contains a propagated discontinuity induced by the dependent delay `lag` and the
discontinuity at time point `T`, or `nothing`.

The interval is estimated by checking the signs of `T + lag(u(t), p, t) - t` for time points
`integrator.t .+ θs` in the interval `[integrator.t, integrator.t + integrator.dt]`.
"""
function discontinuity_interval(integrator::DDEIntegrator, lag, T, Θs)
    # use start and end point of last time interval to check for discontinuities
    previous_condition = discontinuity_function(integrator, lag, T, integrator.t)
    if isapprox(
            previous_condition, 0;
            rtol = integrator.discontinuity_reltol,
            atol = integrator.discontinuity_abstol
        )
        prev_sign = 0
    else
        prev_sign = cmp(previous_condition, zero(previous_condition))
    end
    new_condition = discontinuity_function(integrator, lag, T, integrator.t + integrator.dt)
    new_sign = cmp(new_condition, zero(new_condition))

    # if signs are different we already know that a discontinuity exists
    if prev_sign * new_sign < 0
        return (zero(eltype(Θs)), one(eltype(Θs)))
    end

    # recheck interpolation intervals if no discontinuity found yet
    prev_Θ = zero(eltype(Θs))

    for i in 2:length(Θs)
        # evaluate sign at next interpolation point
        new_Θ = Θs[i]
        new_t = integrator.t + new_Θ * integrator.dt
        new_condition = discontinuity_function(integrator, lag, T, new_t)
        new_sign = cmp(new_condition, zero(new_condition))

        # return estimated interval if we find a root or observe a switch of signs
        if new_sign == 0
            return (new_Θ, new_Θ)
        elseif prev_sign * new_sign < 0
            return (prev_Θ, new_Θ)
        else
            # otherwise update lower estimate of subinterval
            prev_sign = new_sign
            prev_Θ = new_Θ
        end
    end

    return nothing
end

"""
    discontinuity_time(integrator::DDEIntegrator, lag, T, interval)

Estimate time point of the propagated discontinuity induced by the dependent delay
`lag` and the discontinuity at time point `T` inside the `interval` of the current
integration step of the `integrator`.
"""
function discontinuity_time(integrator::DDEIntegrator, lag, T, (bottom_Θ, top_Θ))
    if bottom_Θ == top_Θ
        # in that case we have already found the time point of a discontinuity
        Θ = top_Θ
    else
        # define function for root finding
        zero_func = let integrator = integrator, lag = lag, T = T, t = integrator.t,
                dt = integrator.dt

            (θ, p = nothing) -> discontinuity_function(integrator, lag, T, t + θ * dt)
        end

        Θ = SimpleNonlinearSolve.solve(
            SimpleNonlinearSolve.IntervalNonlinearProblem{false}(
                zero_func,
                (
                    bottom_Θ,
                    top_Θ,
                )
            ),
            SimpleNonlinearSolve.ITP()
        ).left
    end

    return integrator.t + Θ * integrator.dt
end
