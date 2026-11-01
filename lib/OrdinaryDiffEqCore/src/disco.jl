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

_discontinuity_detection_enabled(cache::AbstractControllerCache) =
    cache.controller.basic.discontinuity_detection
_discontinuity_detection_enabled(cache::DummyControllerCache) =
    cache.discontinuity_detection
_discontinuity_detection_enabled(cache::CompositeControllerCache) =
    _discontinuity_detection_enabled(first(cache.caches))

# The interpolation buffer and the condition buffers of a discontinuity problem's zero function
# (`__init` builds them with these, `_disco_zero` derives their types from them).
_disco_u1(u) = u isa AbstractArray ? similar(u) : zero(u)
function _disco_outs(u, cb)
    if cb isa VectorContinuousCallback
        arr = u isa AbstractArray ? similar(u, cb.len) : zeros(typeof(u), cb.len)
        return arr, similar(arr), similar(arr)
    else
        return nothing, zeros(Float64, 1), zeros(Float64, 1)
    end
end

"""
    _disco_zero(disco_prob, integrator, cb)

The zero function of the discontinuity problem of callback `cb`. `__init` stores it behind a
`FunctionWrapper` so that every problem in `disco_probs` has one concrete type whatever the
callback's type, and a `FunctionWrapper` keeps its payload in an untyped field. When the
integrator's callbacks are a tuple they are the objects `__init` built the zero functions with,
so the zero function's type is determined by the integrator and `cb` and is asserted here,
which makes every load and store of its fields in `find_discontinuity` static. Callbacks held
in a vector are wrapped after `__init` (`_despecialize_callbacks!`) and are read untyped.
"""
@inline function _disco_zero(disco_prob, integrator, cb)
    integrator.opts.callback.continuous_callbacks isa Tuple || return disco_prob.f.f.obj.x
    uT = typeof(integrator.u)
    O = Base.promote_op(_disco_outs, uT, typeof(cb))
    U1 = Base.promote_op(_disco_u1, uT)
    # an inferred type that isn't concrete may differ from the stored one: read untyped then
    (isconcretetype(O) && isconcretetype(U1)) || return disco_prob.f.f.obj.x
    Z = zero_func_struct{
        U1, uT, typeof(integrator.dt), typeof(integrator.k),
        typeof(integrator.cache), typeof(integrator.opts.save_idxs),
        typeof(integrator.differential_vars), typeof(cb), fieldtype(O, 1), fieldtype(O, 2),
        typeof(integrator.f), typeof(integrator.tprev), typeof(integrator.p),
    }
    return (disco_prob.f.f.obj::Base.RefValue{Z}).x
end

"""
    _disco_bracket(disco_prob, θ)

`disco_prob` with the bracket `(0, θ)`: the same function and parameters, a tuple `tspan`.
The stored problem's `tspan` is never written, so the bracket is a `Tuple` wherever the
bracketing solver splats it.
"""
@inline function _disco_bracket(disco_prob, θ)
    return IntervalNonlinearProblem{false}(disco_prob.f, (zero(θ), θ), disco_prob.p)
end

function find_discontinuity(integrator)
    cb = integrator.opts.callback
    dt = integrator.dt
    u = integrator.u
    uprev = integrator.uprev
    integrator.is_disco_step = false
    cb === nothing && return one(dt)
    isempty(cb.continuous_callbacks) && return one(dt)
    p = integrator.p
    t = integrator.t
    k = integrator.k
    # the root fraction only places the retried step, so it is found on primal values
    breakpointθ = one(SciMLBase.value(dt))
    disco_probs = get_disco_probs(integrator.controller_cache)
    idx = 1
    addsteps_called = false
    for i in cb.continuous_callbacks
        if (!(i.maybe_discontinuity::Bool))
            continue
        end
        disco_prob = disco_probs[idx]
        disco_zero = _disco_zero(disco_prob, integrator, i)
        _load_disco_state!(disco_zero, dt, uprev, u, k, t, p)
        if (i isa VectorContinuousCallback)
            len_cb = i.len::Int
            i.condition(disco_zero.out_low, uprev, t, integrator)
            i.condition(disco_zero.out_high, u, t + dt, integrator)
            for j in 1:len_cb
                if (disco_zero.out_low[j] * disco_zero.out_high[j] < zero(disco_zero.out_low[j]))
                    if (!addsteps_called)
                        addsteps_called = true
                        _ode_addsteps!(
                            disco_zero.k, disco_zero.tprev, disco_zero.uprev, disco_zero.u,
                            disco_zero.dt, disco_zero.f, disco_zero.p, disco_zero.cache, false, true, false
                        )
                    end
                    _set_disco_index!(disco_zero, j)
                    sol = solve(_disco_bracket(disco_prob, breakpointθ))
                    tmp = sol[]
                    if (sol.retcode != ReturnCode.InitialFailure && !isnan(tmp) && tmp < breakpointθ)
                        breakpointθ = tmp
                        integrator.is_disco_step = true
                        integrator.disco_checkpoint = integrator.t + dt #our prev rejected step, we shouldn't step too far past this
                    end
                end
            end
        else
            disco_zero.out_low[1] = SciMLBase.value(i.condition(uprev, t, integrator))
            disco_zero.out_high[1] = SciMLBase.value(i.condition(u, t + dt, integrator))
            if (disco_zero.out_low[1] * disco_zero.out_high[1] < zero(disco_zero.out_low[1]))
                if (!addsteps_called)
                    addsteps_called = true
                    _ode_addsteps!(
                        disco_zero.k, disco_zero.tprev, disco_zero.uprev, disco_zero.u,
                        disco_zero.dt, disco_zero.f, disco_zero.p, disco_zero.cache, false, true, false
                    )
                end
                sol = solve(_disco_bracket(disco_prob, breakpointθ))
                tmp = sol[]
                if (sol.retcode != ReturnCode.InitialFailure && !isnan(tmp) && tmp < breakpointθ)
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
