abstract type AbstractExpCache{expOpType <: AbstractSciMLOperator} end

struct ExpCache{expOpType} <: AbstractExpCache{expOpType}
    expÂ_cached::Array{expOpType, 2}
    expÂ_for_this_step::Vector{expOpType}
end
struct ExpCacheNoLdiv{expOpType} <: AbstractExpCache{expOpType}
    exp_cache::ExpCache{expOpType}
    exp_cache_inv::ExpCache{expOpType}
end

function get_op_for_this_step(cache::ExpCache{expOpType}, index::Int) where {expOpType}
    cache.expÂ_for_this_step[index]
end
function get_op_for_this_step(cache_no_ldiv::ExpCacheNoLdiv{expOpType},
        positive::Bool, index::Int) where {expOpType}
    positive ? cache_no_ldiv.exp_cache.expÂ_for_this_step[index] :
    cache_no_ldiv.exp_cache_inv.expÂ_for_this_step[index]
end

mutable struct RKIPCache{
    expOpType <: AbstractSciMLOperator, cacheType <: AbstractExpCache{expOpType},
    tType <: Number, opType <: AbstractSciMLOperator, uType, iip} <:
               OrdinaryDiffEqMutableCache
    exp_cache::cacheType
    last_step::tType
    cached::Vector{Bool}
    tmp::uType
    utilde::uType
    kk::Vector{uType}
    c_unique::Vector{tType}
    c_mapping::Vector{Integer}
end

get_fsalfirstlast(cache::RKIPCache, u) = (zero(cache.tmp), zero(cache.tmp))

@inline function cache_exp!(cache::ExpCache{expOpType},
        A::opType,
        h::T,
        action::Symbol,
        step_index::Int,
        unique_stage_index::Int) where {
        expOpType <: AbstractSciMLOperator, opType <: AbstractSciMLOperator, T <: Number}
    (; expÂ_for_this_step, expÂ_cached) = cache
    expÂ_for_this_step[unique_stage_index] = (action == :use_cached) ?
                                              expÂ_cached[step_index, unique_stage_index] :
                                              exp(A, h) # fetching or generating exp(Â*c_i*dt)
    if action == :cache
        expÂ_cached[step_index, unique_stage_index] = expÂ_for_this_step[unique_stage_index] # storing exp(Â*c_i*dt)
    end
end

@inline function cache_exp!(cache::ExpCacheNoLdiv{expOpType},
        Â::opType,
        h::T,
        action::Symbol,
        step_index::Int,
        unique_stage_index::Int) where {
        expOpType <: AbstractSciMLOperator, opType <: AbstractSciMLOperator, T <: Number}
    cache_exp!(cache.exp_cache, Â, h, action, step_index, unique_stage_index)
    cache_exp!(cache.exp_cache_inv, Â, -h, action, step_index, unique_stage_index)
end

"""
    Prepare/generate all the needed exp(± A * dt * c[i]) for a step size dt
"""
@inline function cache_exp_op_for_this_step!(
        cache::RKIPCache{expOpType, cacheType, tType, opType, uType, iip},
        Â::opType, dt::tType,
        alg::algType) where {expOpType, cacheType, tType, opType, uType, algType, iip}
    (; dt_for_expÂ_caching) = alg

    if !iszero(dt) && !(dt ≈ cache.last_step) # we check that new exp(A dt) are needed
        dt_abs = abs(dt) # only the positive dt are used for indexing
        action = :single_use # exp(A*dt) is only computed for this step

        step_index = clamp(searchsortedlast(dt_for_expÂ_caching, dt_abs),
            1, lastindex(dt_for_expÂ_caching)) # fetching the index corresponding to the step size

        if dt_for_expÂ_caching[step_index] ≈ dt_abs # if dt corresponds to a cahing step
            action = (cache.cached[step_index] ? :use_cached : :cache) # if alreay present, we reuse the cached, otherwise it is generated
        end

        for (unique_stage_index, c) in enumerate(cache.c_unique) # iterating over all unique c_i of the RK tableau
            cache_exp!(
                cache.exp_cache, Â, abs(dt * c), action, step_index, unique_stage_index) # generating and caching
        end
        cache.cached[step_index] |= (action == :cache) # set the flag that we have already cached exp(Â*c_i*dt) for this dt
    end

    cache.last_step = dt
end
