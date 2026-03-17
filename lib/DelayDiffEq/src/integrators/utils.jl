"""
    advance_or_update_ode_integrator!(integrator::DDEIntegrator[, always_calc_begin = false])

Advance or update the ODE integrator of `integrator` to the next time interval by updating its
values and interpolation data with the current values and a full set of interpolation data of
`integrator`.
"""
function advance_or_update_ode_integrator!(integrator, always_calc_begin = false)
    ode_integrator = integrator.integrator

    return if ode_integrator.t != integrator.t + integrator.dt
        advance_ode_integrator!(integrator, always_calc_begin)
    else
        update_ode_integrator!(integrator, always_calc_begin)
    end
end

# Helper function to get the current cache from composite algorithm caches
# Handles both regular cache.caches arrays and DefaultCache with individual fields
function get_current_cache(cache, current)
    if cache isa OrdinaryDiffEqCore.DefaultCache
        # DefaultCache stores caches as individual fields (cache1, cache2, etc.)
        # Use getfield to safely access the cache fields
        if current == 1
            return cache.cache1
        elseif current == 2
            return cache.cache2
        elseif current == 3
            return isdefined(cache, :cache3) ? cache.cache3 : cache.cache1
        elseif current == 4
            return isdefined(cache, :cache4) ? cache.cache4 : cache.cache1
        elseif current == 5
            return isdefined(cache, :cache5) ? cache.cache5 : cache.cache1
        else
            return isdefined(cache, :cache6) ? cache.cache6 : cache.cache1
        end
    else
        return cache.caches[current]
    end
end

"""
    advance_ode_integrator!(integrator::DDEIntegrator[, always_calc_begin = false])

Advance the ODE integrator of `integrator` to the next time interval by updating its values
and interpolation data with the current values and a full set of interpolation data of
`integrator`.
"""
function advance_ode_integrator!(integrator::DDEIntegrator, always_calc_begin = false)
    (; f, u, t, p, k, dt, uprev, alg, cache) = integrator
    ode_integrator = integrator.integrator

    # algorithm only works if current time of DDE integrator equals final time point
    # of solution
    t != ode_integrator.sol.t[end] && error("cannot advance ODE integrator")

    # complete interpolation data of DDE integrator for time interval [t, t+dt]
    # and copy it to ODE integrator
    # has to be done before updates to ODE integrator, otherwise history function
    # is incorrect
    if iscomposite(alg)
        OrdinaryDiffEqCore._ode_addsteps!(
            k, t, uprev, u, dt, f, p, get_current_cache(cache, cache.current),
            always_calc_begin, true, true
        )
    else
        OrdinaryDiffEqCore._ode_addsteps!(
            k, t, uprev, u, dt, f, p, cache, always_calc_begin,
            true, true
        )
    end
    @inbounds for i in 1:length(k)
        copyat_or_push!(ode_integrator.k, i, k[i])
    end

    # move ODE integrator to interval [t, t+dt]
    ode_integrator.t = t + dt
    ode_integrator.tprev = t
    ode_integrator.dt = dt
    if iscomposite(alg)
        ode_integrator.cache.current = cache.current
    end
    if isinplace(integrator.sol.prob)
        recursivecopy!(ode_integrator.u, u)
    else
        ode_integrator.u = u
    end

    # u(t) is not modified hence we do not have to copy it
    ode_integrator.uprev = ode_integrator.sol.u[end]

    # update prev_idx to index of t and u(t) in solution
    integrator.prev_idx = length(ode_integrator.sol.t)

    return nothing
end

"""
    update_ode_integrator!(integrator::DDEIntegrator[, always_calc_begin = false])

Update the ODE integrator of `integrator` by updating its values and interpolation data
with the current values and a full set of interpolation data of `integrator`.
"""
function update_ode_integrator!(integrator::DDEIntegrator, always_calc_begin = false)
    (; f, u, t, p, k, dt, uprev, alg, cache) = integrator
    ode_integrator = integrator.integrator

    # algorithm only works if the ODE integrator is already moved to the current integration
    # interval
    ode_integrator.t != t + dt && error("cannot update ODE integrator")

    if iscomposite(alg)
        OrdinaryDiffEqCore._ode_addsteps!(
            k, t, uprev, u, dt, f, p, get_current_cache(cache, cache.current),
            always_calc_begin, true, true
        )
    else
        OrdinaryDiffEqCore._ode_addsteps!(
            k, t, uprev, u, dt, f, p, cache,
            always_calc_begin, true, true
        )
    end
    @inbounds for i in 1:length(k)
        copyat_or_push!(ode_integrator.k, i, k[i])
    end

    # update state of the dummy ODE solver
    if isinplace(integrator.sol.prob)
        recursivecopy!(ode_integrator.u, u)
    else
        ode_integrator.u = integrator.u
    end

    return nothing
end

"""
    move_back_ode_integrator!(integrator::DDEIntegrator)

Move the ODE integrator of `integrator` one integration step back by reverting its values
and interpolation data to the values saved in the dense history.
"""
function move_back_ode_integrator!(integrator::DDEIntegrator)
    ode_integrator = integrator.integrator
    (; sol) = ode_integrator

    # set values of the ODE integrator back to the values in the solution
    if isinplace(sol.prob)
        recursivecopy!(ode_integrator.u, sol.u[end])
    else
        ode_integrator.u = sol.u[end]
    end
    ode_integrator.t = sol.t[end]
    ode_integrator.tprev = sol.t[integrator.prev2_idx]

    # u(tprev) is not modified hence we do not have to copy it
    ode_integrator.uprev = sol.u[integrator.prev2_idx]

    # revert to the previous time step
    ode_integrator.dt = ode_integrator.dtcache

    # we do not have to reset the interpolation data in the initial time step since always a
    # constant extrapolation is used (and interpolation data of solution at initial
    # time point is not complete!)
    if length(sol.t) > 1
        recursivecopy!(ode_integrator.k, sol.k[end])
    end

    return nothing
end

#=
Dealing with discontinuities

If we hit a discontinuity (this is checked in `apply_step!`), then we remove the
discontinuity, additional discontinuities at the current time point (if present), and
maybe also discontinuities and time stops coming shortly after the current time point
in `handle_discontinuities!`. The order of the discontinuity at the current time point is
defined as the lowest order of all these discontinuities.

If the problem is not neutral, we will only add additional discontinuities if
this order is less or equal to the order of the algorithm in
`add_next_discontinuities!`. If we add discontinuities, we add discontinuities
of the next order caused by constant lags (these we can calculate explicitly and
just add them to `d_discontinuities` and `tstops`) and we add the current
discontinuity to `tracked_discontinuities` which is the array of old
discontinuities that are checked by a `DiscontinuityCallback` (if existent).
=#

# handle discontinuities at the current time point of the `integrator`
function OrdinaryDiffEqCore.handle_discontinuities!(integrator::DDEIntegrator)
    # remove all discontinuities at current time point and calculate minimal order
    # of these discontinuities
    d = OrdinaryDiffEqCore.pop_discontinuity!(integrator)
    order = d.order
    tdir_t = integrator.tdir * integrator.t

    @SciMLMessage(
        lazy"Handling discontinuity at t = $(integrator.t) with order $order",
        integrator.opts.verbose, :discontinuity_tracking
    )

    while OrdinaryDiffEqCore.has_discontinuity(integrator) &&
            OrdinaryDiffEqCore.first_discontinuity(integrator) == tdir_t
        d2 = OrdinaryDiffEqCore.pop_discontinuity!(integrator)
        order = min(order, d2.order)
    end

    # remove all discontinuities close to the current time point as well and
    # calculate minimal order of these discontinuities
    # integrator.EEst has unitless type of integrator.t
    if integrator.EEst isa AbstractFloat
        maxΔt = 10eps(integrator.t)

        while OrdinaryDiffEqCore.has_discontinuity(integrator) &&
                abs(OrdinaryDiffEqCore.first_discontinuity(integrator).t - tdir_t) < maxΔt
            d2 = OrdinaryDiffEqCore.pop_discontinuity!(integrator)
            order = min(order, d2.order)
        end

        # also remove all corresponding time stops
        while OrdinaryDiffEqCore.has_tstop(integrator) &&
                abs(OrdinaryDiffEqCore.first_tstop(integrator) - tdir_t) < maxΔt
            OrdinaryDiffEqCore.pop_tstop!(integrator)
        end
    end

    # add discontinuities of next order to integrator
    add_next_discontinuities!(integrator, order)

    return nothing
end

"""
    add_next_discontinuities!(integrator::DDEIntegrator, order[, t=integrator.t])

Add discontinuities of next order that are propagated from discontinuity of
order `order` at time `t` in `integrator`, but only if `order` is less or equal
than the order of the applied method or the problem is neutral.

Discontinuities caused by constant delays are immediately calculated, and
discontinuities caused by dependent delays are tracked by a callback.
"""
function add_next_discontinuities!(integrator, order, t = integrator.t)
    neutral = integrator.sol.prob.neutral
    next_order = neutral ? order : order + 1

    if neutral
        @SciMLMessage(
            lazy"Adding next discontinuities for neutral DDE: order remains $order",
            integrator.opts.verbose, :neutral_delay
        )
    else
        @SciMLMessage(
            lazy"Adding next discontinuities: order increases from $order to $next_order",
            integrator.opts.verbose, :discontinuity_tracking
        )
    end

    # only track discontinuities up to order of the applied method
    alg_maximum_order = OrdinaryDiffEqCore.alg_maximum_order(integrator.alg)
    next_order <= alg_maximum_order + 1 || return

    # discontinuities caused by constant lags
    if has_constant_lags(integrator)
        constant_lags = integrator.sol.prob.constant_lags
        maxlag = integrator.tdir * (integrator.sol.prob.tspan[end] - t)

        for lag in constant_lags
            if integrator.tdir * lag < maxlag
                # calculate discontinuity and add it to heap of discontinuities and time stops
                d = Discontinuity(integrator.tdir * (t + lag), next_order)
                push!(integrator.d_discontinuities_propagated, d)
                push!(integrator.tstops_propagated, d.t)
            end
        end
    end

    # track propagated discontinuities with callback
    push!(integrator.tracked_discontinuities, Discontinuity(integrator.tdir * t, order))

    return nothing
end

# Interface for accessing and removing next time stops and discontinuities
function OrdinaryDiffEqCore.has_tstop(integrator::DDEIntegrator)
    return _has(integrator.opts.tstops, integrator.tstops_propagated)
end
function OrdinaryDiffEqCore.first_tstop(integrator::DDEIntegrator)
    return _first(integrator.opts.tstops, integrator.tstops_propagated)
end
function OrdinaryDiffEqCore.pop_tstop!(integrator::DDEIntegrator)
    return _pop!(integrator.opts.tstops, integrator.tstops_propagated)
end

function OrdinaryDiffEqCore.has_discontinuity(integrator::DDEIntegrator)
    return _has(
        integrator.opts.d_discontinuities,
        integrator.d_discontinuities_propagated
    )
end
function OrdinaryDiffEqCore.first_discontinuity(integrator::DDEIntegrator)
    return _first(
        integrator.opts.d_discontinuities,
        integrator.d_discontinuities_propagated
    )
end
function OrdinaryDiffEqCore.pop_discontinuity!(integrator::DDEIntegrator)
    return _pop!(
        integrator.opts.d_discontinuities,
        integrator.d_discontinuities_propagated
    )
end

_has(x, y) = !isempty(x) || !isempty(y)
function _first(x, y)
    if isempty(x)
        return first(y)
    elseif isempty(y)
        return first(x)
    else
        return min(first(x), first(y))
    end
end
function _pop!(x, y)
    return if isempty(x)
        pop!(y)
    elseif isempty(y)
        pop!(x)
    else
        if first(x) < first(y)
            pop!(x)
        else
            pop!(y)
        end
    end
end
