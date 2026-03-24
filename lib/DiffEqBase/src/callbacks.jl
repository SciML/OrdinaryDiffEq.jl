"""
    initialize!(cb::CallbackSet,u,t,integrator::DEIntegrator)

Recursively apply `initialize!` and return whether any modified u
"""
function initialize!(cb::CallbackSet, u, t, integrator::DEIntegrator)
    return initialize!(
        u, t, integrator, false, cb.continuous_callbacks...,
        cb.discrete_callbacks...
    )
end
initialize!(cb::CallbackSet{Tuple{}, Tuple{}}, u, t, integrator::DEIntegrator) = false
function initialize!(
        u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback, cs::DECallback...
    )
    c.initialize(c, u, t, integrator)
    return initialize!(u, t, integrator, any_modified || integrator.u_modified, cs...)
end
function initialize!(
        u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback
    )
    c.initialize(c, u, t, integrator)
    return any_modified || integrator.u_modified
end

"""
    finalize!(cb::CallbackSet,u,t,integrator::DEIntegrator)

Recursively apply `finalize!` and return whether any modified u
"""
function finalize!(cb::CallbackSet, u, t, integrator::DEIntegrator)
    return finalize!(u, t, integrator, false, cb.continuous_callbacks..., cb.discrete_callbacks...)
end
finalize!(cb::CallbackSet{Tuple{}, Tuple{}}, u, t, integrator::DEIntegrator) = false
function finalize!(
        u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback, cs::DECallback...
    )
    c.finalize(c, u, t, integrator)
    return finalize!(u, t, integrator, any_modified || integrator.u_modified, cs...)
end
function finalize!(
        u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback
    )
    c.finalize(c, u, t, integrator)
    return any_modified || integrator.u_modified
end

# Helpers
function Base.isempty(cb::CallbackSet)
    return isempty(cb.continuous_callbacks) && isempty(cb.discrete_callbacks)
end
Base.isempty(cb::AbstractContinuousCallback) = false
Base.isempty(cb::AbstractDiscreteCallback) = false

has_continuous_callback(cb::DiscreteCallback) = false
has_continuous_callback(cb::ContinuousCallback) = true
has_continuous_callback(cb::VectorContinuousCallback) = true
has_continuous_callback(cb::CallbackSet) = !isempty(cb.continuous_callbacks)
has_continuous_callback(cb::Nothing) = false

rightfloat(t, tdir) = isone(tdir) ? nextfloat(t) : prevfloat(t)

# Callback handling

function get_tmp(integrator::DEIntegrator, callback)
    _tmp = get_tmp_cache(integrator)
    _tmp === nothing && return nothing
    _cache = first(_tmp)
    if callback.idxs === nothing
        tmp = _cache
    elseif !(callback.idxs isa Number)
        tmp = @view _cache[callback.idxs]
    else
        tmp = nothing
    end
    return tmp
end

function get_condition(integrator::DEIntegrator, callback, abst)
    tmp = get_tmp(integrator, callback)
    ismutable = !(tmp === nothing)
    if abst == integrator.t
        if callback.idxs === nothing
            tmp = integrator.u
        elseif callback.idxs isa Number
            tmp = integrator.u[callback.idxs]
        else
            tmp = @view integrator.u[callback.idxs]
        end
    elseif abst == integrator.tprev
        if callback.idxs === nothing
            tmp = integrator.uprev
        elseif callback.idxs isa Number
            tmp = integrator.uprev[callback.idxs]
        else
            tmp = @view integrator.uprev[callback.idxs]
        end
    else
        if ismutable
            if callback.idxs === nothing
                integrator(tmp, abst, Val{0})
            else
                integrator(tmp, abst, Val{0}, idxs = callback.idxs)
            end
        else
            if callback.idxs === nothing
                tmp = integrator(abst, Val{0})
            else
                tmp = integrator(abst, Val{0}, idxs = callback.idxs)
            end
        end
        # ismutable && !(callback.idxs isa Number) ? integrator(tmp,abst,Val{0},idxs=callback.idxs) :
        #                                                 tmp = integrator(abst,Val{0},idxs=callback.idxs)
    end
    integrator.sol.stats.ncondition += 1
    if callback isa VectorContinuousCallback
        callback.condition(
            @view(integrator.callback_cache.tmp_condition[1:(callback.len)]),
            tmp, abst, integrator
        )
        return @view(integrator.callback_cache.tmp_condition[1:(callback.len)])
    else
        return callback.condition(tmp, abst, integrator)
    end
end

# Use a generated function for type stability even when many callbacks are given
@inline function find_first_continuous_callback(
        integrator,
        callbacks::Vararg{
            AbstractContinuousCallback,
            N,
        }
    ) where {N}
    return find_first_continuous_callback(integrator, tuple(callbacks...))
end
@generated function find_first_continuous_callback(
        integrator,
        callbacks::NTuple{
            N,
            AbstractContinuousCallback,
        }
    ) where {N}
    ex = quote
        tmin, upcrossing,
            event_occurred, event_idx, residual = find_callback_time(
            integrator,
            callbacks[1], 1
        )
        identified_idx = 1
    end
    for i in 2:N
        ex = quote
            $ex
            tmin2, upcrossing2,
                event_occurred2, event_idx2, residual2 = find_callback_time(
                integrator,
                callbacks[$i],
                $i
            )
            if event_occurred2 && (!event_occurred || integrator.tdir * tmin2 < integrator.tdir * tmin)
                tmin = tmin2
                upcrossing = upcrossing2
                event_occurred = true
                event_idx = event_idx2
                identified_idx = $i
                residual = residual2
            end
        end
    end
    ex = quote
        $ex
        if event_occurred
            integrator.last_event_error = value(residual)
        end
        return tmin, upcrossing, event_occurred, event_idx, identified_idx, $N
    end
    return ex
end

@inline function find_callback_time(
        integrator, callback::VectorContinuousCallback,
        callback_idx
    )
    if callback.interp_points != 0
        addsteps!(integrator)
    end

    # Compute previous sign
    bottom_sign = @view(integrator.callback_cache.prev_sign[1:(callback.len)])
    bottom_t = integrator.tprev
    bottom_condition = get_condition(integrator, callback, integrator.tprev)
    @. bottom_sign = sign(bottom_condition)

    if integrator.event_last_time == callback_idx
        nudged_idx = integrator.vector_event_last_time
        # If there was a previous event, nudge tprev on the right
        # side of the root (if necessary) to avoid repeat detection

        if callback.interp_points == 0
            addsteps!(integrator)
        end

        # Evaluate condition slightly in future
        nudged_t = nudge_tprev(integrator, callback, ArrayInterface.allowed_getindex(bottom_condition, nudged_idx))
        tmp_condition = get_condition(integrator, callback, nudged_t)

        ArrayInterface.allowed_setindex!(bottom_sign, sign(ArrayInterface.allowed_getindex(tmp_condition, nudged_idx)), nudged_idx)
    else
        nudged_idx = -1
        nudged_t = bottom_t
    end

    # Check if an event occured
    event_occurred, event_idx, top_t, top_sign =
        check_event_occurence(integrator, callback, bottom_sign)

    # Find callback time if occurence
    if !event_occurred
        callback_t = integrator.t
        min_event_idx = 1
        residual = zero(eltype(bottom_condition))
    elseif isdiscrete(integrator.alg) || callback.rootfind == SciMLBase.NoRootFind
        callback_t = top_t
        min_event_idx = 1
        for i in 1:length(event_idx)
            if ArrayInterface.allowed_getindex(event_idx, i) == 1
                min_event_idx = i
                break
            end
        end
        residual = zero(eltype(bottom_condition))
    else
        callback_t = rightfloat(top_t, integrator.tdir)
        min_event_idx = -1
        for idx in 1:length(event_idx)
            if ArrayInterface.allowed_getindex(event_idx, idx) != 0
                function zero_func(abst, p = nothing)
                    return ArrayInterface.allowed_getindex(
                        get_condition(
                            integrator,
                            callback,
                            abst
                        ), idx
                    )
                end
                if iszero(ArrayInterface.allowed_getindex(top_sign, idx))
                    cbi_t = top_t
                else
                    if idx == nudged_idx
                        cbi_t = find_root(zero_func, (nudged_t, top_t), callback.rootfind)
                    else
                        cbi_t = find_root(zero_func, (bottom_t, top_t), callback.rootfind)
                    end
                end
                if integrator.tdir * cbi_t < integrator.tdir * callback_t
                    min_event_idx = idx
                    callback_t = cbi_t
                    residual = zero_func(cbi_t)
                end
            end
        end

        if min_event_idx < 0
            error("Callback handling failed. Please file an issue with code to reproduce.")
        end
    end

    return callback_t, ArrayInterface.allowed_getindex(bottom_sign, min_event_idx),
        event_occurred::Bool, min_event_idx::Int, residual
end

@inline function find_callback_time(
        integrator, callback::ContinuousCallback,
        callback_idx
    )
    if callback.interp_points != 0
        addsteps!(integrator)
    end

    # Compute previous sign
    bottom_t = integrator.tprev
    bottom_condition = get_condition(integrator, callback, bottom_t)
    if integrator.event_last_time == callback_idx
        # If there was a previous event, nudge tprev on the right
        # side of the root (if necessary) to avoid repeat detection

        if callback.interp_points == 0
            addsteps!(integrator)
        end

        bottom_t = nudge_tprev(integrator, callback, bottom_condition)
        bottom_condition = get_condition(integrator, callback, bottom_t)
    end
    bottom_sign = sign(bottom_condition)

    # Check if an event occured
    event_occurred, event_idx, top_t, top_sign =
        check_event_occurence(integrator, callback, bottom_sign)

    if !event_occurred
        callback_t = integrator.t
        residual = zero(bottom_condition)
    elseif isdiscrete(integrator.alg) || callback.rootfind == SciMLBase.NoRootFind || iszero(top_sign)
        callback_t = top_t
        residual = zero(bottom_condition)
    else
        # Find callback time
        zero_func(abst, p = nothing) = get_condition(integrator, callback, abst)
        callback_t = find_root(zero_func, (bottom_t, top_t), callback.rootfind)
        residual = zero_func(callback_t)
    end

    return callback_t, bottom_sign, event_occurred, event_idx, residual
end

"""
Return a nudged (if neccessary) value of `integrator.tprev` to avoid repeat event detection
- `integrator`
- `callback`: Last occuring callback
- `condition_tprev`: Condition of last occuring callback evaluated at `integrator.tprev`
"""
function nudge_tprev(integrator, callback, condition_tprev)
    # Assume the previous event might affect the condition/root
    return if abs(condition_tprev - integrator.last_event_error) <= callback.abstol
        # We are still close to the root
        right_t = integrator.tprev + integrator.dt * callback.repeat_nudge
    else
        # We are far away from the root, keep the current sign
        right_t = integrator.tprev
    end
end

"""
Determine if an event occured in the integration time step
"""
function check_event_occurence(integrator, callback, bottom_sign)
    top_t = integrator.t
    event_occurred, event_idx, top_sign =
        check_event_occurence_upto(integrator, callback, bottom_sign, top_t)

    if callback.interp_points != 0 && !isdiscrete(integrator.alg) &&
            any(iszero, event_idx)
        # Use the interpolants for safety checking
        ts = range(integrator.tprev, stop = integrator.t, length = callback.interp_points)
        for i in 2:length(ts)
            top_t = ts[i]
            event_occurred, event_idx, top_sign =
                check_event_occurence_upto(integrator, callback, bottom_sign, top_t)
            if event_occurred
                break
            end
        end
    end

    return event_occurred, event_idx, top_t, top_sign
end

"""
Determine if an event occured before `top_t``
"""
function check_event_occurence_upto(integrator, callback::ContinuousCallback, bottom_sign, top_t)
    top_sign = sign(get_condition(integrator, callback, top_t))
    event_occurred = is_event_occurence(bottom_sign, top_sign, callback.affect!, callback.affect_neg!)
    event_idx = event_occurred ? 1.0 : 0.0
    return event_occurred, event_idx, top_sign
end

function check_event_occurence_upto(integrator, callback::VectorContinuousCallback, bottom_sign, top_t)
    event_idx = top_condition = @views(integrator.callback_cache.next_condition[1:(callback.len)])
    top_sign = @view(integrator.callback_cache.next_sign[1:(callback.len)])
    copyto!(top_condition, get_condition(integrator, callback, top_t))
    @. top_sign = sign(top_condition)

    # Determine event occurence
    event_occurred = findall_events!(
        top_condition, callback.affect!, callback.affect_neg!,
        bottom_sign
    )
    return event_occurred, event_idx, top_sign
end

"""
Find either exact or floating point precision root of `f`.
If the exact root cannot be represented, return closest floating point number depending on `rootfind`

Assumes that:
 - `tup[1] < tup[2]` for a forward integration
 - `tup[1] > tup[2]` for a backward integration
 - The nonlinear solver return left/right roots in the same order as `tup[1]`/`tup[2]`
"""
function find_root(f, tup, rootfind::SciMLBase.RootfindOpt)
    sol = solve(
        IntervalNonlinearProblem{false}(f, tup),
        ModAB(), abstol = 0.0, reltol = 0.0
    )
    if rootfind == SciMLBase.LeftRootFind
        return sol.left
    else
        return sol.right
    end
end

"""
findall_events!(next_sign,affect!,affect_neg!,prev_sign)

Modifies `next_sign` to be an array of booleans for if there is a sign change
in the interval between prev_sign and next_sign.
Return `true` if any event occured.
"""
function findall_events!(
        next_sign::Union{Array, SubArray}, affect!::F1, affect_neg!::F2,
        prev_sign::Union{Array, SubArray}
    ) where {F1, F2}
    @inbounds for i in 1:length(prev_sign)
        next_sign[i] = (
            (prev_sign[i] < 0 && affect! !== nothing) ||
                (prev_sign[i] > 0 && affect_neg! !== nothing)
        ) &&
            prev_sign[i] * next_sign[i] <= 0
    end
    return any(isone, next_sign)
end

function findall_events!(next_sign, affect!::F1, affect_neg!::F2, prev_sign) where {F1, F2}
    hasaffect::Bool = affect! !== nothing
    hasaffectneg::Bool = affect_neg! !== nothing
    f = (n, p) -> ((p < 0 && hasaffect) || (p > 0 && hasaffectneg)) && p * n <= 0
    A = map!(f, next_sign, next_sign, prev_sign)
    return any(isone, next_sign)
end

"""
Return `true` if an event occured.
"""
function is_event_occurence(prev_sign::Number, next_sign::Number, affect!::F1, affect_neg!::F2) where {F1, F2}
    return (
        (prev_sign < 0 && affect! !== nothing) ||
            (prev_sign > 0 && affect_neg! !== nothing)
    ) && prev_sign * next_sign <= 0
end

function apply_callback!(
        integrator,
        callback::Union{ContinuousCallback, VectorContinuousCallback},
        cb_time, prev_sign, event_idx
    )
    if isadaptive(integrator)
        set_proposed_dt!(
            integrator,
            integrator.tdir * max(
                nextfloat(integrator.opts.dtmin),
                integrator.tdir * callback.dtrelax * integrator.dt
            )
        )
    end

    change_t_via_interpolation!(
        integrator, cb_time, Val{:false}, callback.initializealg
    )

    # handle saveat
    _, savedexactly = savevalues!(integrator)
    saved_in_cb = true

    @inbounds if callback.save_positions[1]
        # if already saved then skip saving
        savedexactly || savevalues!(integrator, true)
    end

    integrator.u_modified = true

    if prev_sign < 0
        if callback.affect! === nothing
            integrator.u_modified = false
        else
            callback isa VectorContinuousCallback ?
                callback.affect!(integrator, event_idx) : callback.affect!(integrator)
        end
    elseif prev_sign > 0
        if callback.affect_neg! === nothing
            integrator.u_modified = false
        else
            callback isa VectorContinuousCallback ?
                callback.affect_neg!(integrator, event_idx) : callback.affect_neg!(integrator)
        end
    end

    if integrator.u_modified
        reeval_internals_due_to_modification!(
            integrator, callback_initializealg = callback.initializealg
        )

        @inbounds if callback.save_positions[2]
            savevalues!(integrator, true)
            if !isdefined(integrator.opts, :save_discretes) || integrator.opts.save_discretes
                if callback isa VectorContinuousCallback
                    SciMLBase.save_discretes!(integrator, callback, event_idx)
                else
                    SciMLBase.save_discretes!(integrator, callback)
                end
            end
            saved_in_cb = true
        end
        return true, saved_in_cb
    end
    return false, saved_in_cb
end

#Base Case: Just one
@inline function apply_discrete_callback!(integrator, callback::DiscreteCallback)
    saved_in_cb = false
    if callback.condition(integrator.u, integrator.t, integrator)
        # handle saveat
        _, savedexactly = savevalues!(integrator)
        saved_in_cb = true
        @inbounds if callback.save_positions[1]
            # if already saved then skip saving
            savedexactly || savevalues!(integrator, true)
        end
        integrator.u_modified = true
        callback.affect!(integrator)
        if integrator.u_modified
            reeval_internals_due_to_modification!(
                integrator, false, callback_initializealg = callback.initializealg
            )
        end
        @inbounds if callback.save_positions[2]
            savevalues!(integrator, true)
            if !isdefined(integrator.opts, :save_discretes) || integrator.opts.save_discretes
                SciMLBase.save_discretes!(integrator, callback)
            end
            saved_in_cb = true
        end
    end
    integrator.sol.stats.ncondition += 1
    return integrator.u_modified, saved_in_cb
end

#Starting: Get bool from first and do next
@inline function apply_discrete_callback!(integrator, callback::DiscreteCallback, args...)
    return apply_discrete_callback!(
        integrator, apply_discrete_callback!(integrator, callback)...,
        args...
    )
end

@inline function apply_discrete_callback!(
        integrator, discrete_modified::Bool,
        saved_in_cb::Bool, callback::DiscreteCallback,
        args...
    )
    bool,
        saved_in_cb2 = apply_discrete_callback!(
        integrator,
        apply_discrete_callback!(
            integrator,
            callback
        )...,
        args...
    )
    return discrete_modified || bool, saved_in_cb || saved_in_cb2
end

@inline function apply_discrete_callback!(
        integrator, discrete_modified::Bool,
        saved_in_cb::Bool, callback::DiscreteCallback
    )
    bool, saved_in_cb2 = apply_discrete_callback!(integrator, callback)
    return discrete_modified || bool, saved_in_cb || saved_in_cb2
end

function max_vector_callback_length_int(cs::CallbackSet)
    return max_vector_callback_length_int(cs.continuous_callbacks...)
end
max_vector_callback_length_int() = nothing
function max_vector_callback_length_int(continuous_callbacks...)
    all(cb -> cb isa ContinuousCallback, continuous_callbacks) && return nothing
    maxlen = -1
    for cb in continuous_callbacks
        if cb isa VectorContinuousCallback && cb.len > maxlen
            maxlen = cb.len
        end
    end
    return maxlen
end

function max_vector_callback_length(cs::CallbackSet)
    continuous_callbacks = cs.continuous_callbacks
    maxlen_cb = nothing
    maxlen = -1
    for cb in continuous_callbacks
        if cb isa VectorContinuousCallback && cb.len > maxlen
            maxlen = cb.len
            maxlen_cb = cb
        end
    end
    return maxlen_cb
end

"""
$(TYPEDEF)
"""
mutable struct CallbackCache{conditionType, signType}
    tmp_condition::conditionType
    next_condition::conditionType
    next_sign::signType
    prev_sign::signType
end

function CallbackCache(
        u, max_len, ::Type{conditionType},
        ::Type{signType}
    ) where {conditionType, signType}
    tmp_condition = similar(u, conditionType, max_len)
    next_condition = similar(u, conditionType, max_len)
    next_sign = similar(u, signType, max_len)
    prev_sign = similar(u, signType, max_len)
    return CallbackCache(tmp_condition, next_condition, next_sign, prev_sign)
end

function CallbackCache(
        max_len, ::Type{conditionType},
        ::Type{signType}
    ) where {conditionType, signType}
    tmp_condition = zeros(conditionType, max_len)
    next_condition = zeros(conditionType, max_len)
    next_sign = zeros(signType, max_len)
    prev_sign = zeros(signType, max_len)
    return CallbackCache(tmp_condition, next_condition, next_sign, prev_sign)
end
