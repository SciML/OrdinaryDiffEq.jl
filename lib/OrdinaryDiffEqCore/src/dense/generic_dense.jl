const DERIVATIVE_ORDER_NOT_POSSIBLE_MESSAGE = """
Derivative order too high for interpolation order. An interpolation derivative is
only accurate to a certain derivative. For example, a second order interpolation
is a quadratic polynomial, and thus third derivatives cannot be computed (will be
incorrectly zero). Thus use a solver with a higher order interpolation or compute
the higher order derivative through other means.

You can find the list of available ODE/DAE solvers with their documented interpolations at:

* https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/
* https://docs.sciml.ai/DiffEqDocs/stable/solvers/dae_solve/
"""

using SciMLBase: SENSITIVITY_INTERP_MESSAGE

struct DerivativeOrderNotPossibleError <: Exception end

function Base.showerror(io::IO, e::DerivativeOrderNotPossibleError)
    print(io, DERIVATIVE_ORDER_NOT_POSSIBLE_MESSAGE)
    return println(io, VERBOSE_MSG)
end

## Integrator Dispatches

# Can get rid of an allocation here with a function
# get_tmp_arr(integrator.cache) which gives a pointer to some
# cache array which can be modified.

@inline function _searchsortedfirst(v::AbstractVector, x, lo::Integer, forward::Bool)
    u = oftype(lo, 1)
    lo = lo - u
    hi = length(v) + u
    @inbounds while lo < hi - u
        m = (lo + hi) >>> 1
        @inbounds if (forward && v[m] < x) || (!forward && v[m] > x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

@inline function _searchsortedlast(v::AbstractVector, x, lo::Integer, forward::Bool)
    u = oftype(lo, 1)
    lo = lo - u
    hi = length(v) + u
    @inbounds while lo < hi - u
        m = (lo + hi) >>> 1
        @inbounds if (forward && v[m] > x) || (!forward && v[m] < x)
            hi = m
        else
            lo = m
        end
    end
    return lo
end

@inline function ode_addsteps!(
        integrator, f = integrator.f, always_calc_begin = false,
        allow_calc_end = true, force_calc_end = false
    )
    cache = integrator.cache
    if cache isa CompositeCache
        cache_current = cache.current
        if cache_current == 1
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.caches[1],
                always_calc_begin, allow_calc_end, force_calc_end
            )
        elseif cache_current == 2
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.caches[2],
                always_calc_begin, allow_calc_end, force_calc_end
            )
        else
            @assert length(integrator.cache.caches) >= cache_current
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.caches[cache_current],
                always_calc_begin, allow_calc_end, force_calc_end
            )
        end
    elseif cache isa DefaultCache
        cache_current = cache.current
        if cache_current == 1
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.cache1,
                always_calc_begin, allow_calc_end, force_calc_end
            )
        elseif cache_current == 2
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.cache2,
                always_calc_begin, allow_calc_end, force_calc_end
            )
        elseif cache_current == 3
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.cache3,
                always_calc_begin, allow_calc_end, force_calc_end
            )
        elseif cache_current == 4
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.cache4,
                always_calc_begin, allow_calc_end, force_calc_end
            )
        elseif cache_current == 5
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.cache5,
                always_calc_begin, allow_calc_end, force_calc_end
            )
        elseif cache_current == 6
            _ode_addsteps!(
                integrator.k, integrator.tprev, integrator.uprev,
                integrator.u,
                integrator.dt, f, integrator.p,
                cache.cache6,
                always_calc_begin, allow_calc_end, force_calc_end
            )
        end
    else
        _ode_addsteps!(
            integrator.k, integrator.tprev, integrator.uprev, integrator.u,
            integrator.dt, f, integrator.p, cache,
            always_calc_begin, allow_calc_end, force_calc_end
        )
    end
    return nothing
end
@inline function SciMLBase.addsteps!(integrator::ODEIntegrator, args...)
    return ode_addsteps!(integrator, args...)
end

@inline function ode_interpolant(Θ, integrator::SciMLBase.DEIntegrator, idxs, deriv)
    SciMLBase.addsteps!(integrator)
    if integrator.cache isa CompositeCache
        val = composite_ode_interpolant(
            Θ, integrator, integrator.cache.caches,
            integrator.cache.current, idxs, deriv
        )
    elseif integrator.cache isa DefaultCache
        val = default_ode_interpolant(
            Θ, integrator, integrator.cache,
            integrator.cache.current, idxs, deriv
        )
    else
        val = ode_interpolant(
            Θ, integrator.dt, integrator.uprev, integrator.u,
            integrator.k, integrator.cache, idxs, deriv, integrator.differential_vars
        )
    end
    return val
end

function default_ode_interpolant(
        Θ, integrator, cache::DefaultCache, alg_choice, idxs, deriv
    )
    if alg_choice == 1
        return ode_interpolant(
            Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache1, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 2
        return ode_interpolant(
            Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache2, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 3
        return ode_interpolant(
            Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache3, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 4
        return ode_interpolant(
            Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache4, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 5
        return ode_interpolant(
            Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache5, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 6
        return ode_interpolant(
            Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache6, idxs,
            deriv, integrator.differential_vars
        )
    else
        error("DefaultCache invalid alg_choice. File an issue.")
    end
end

@generated function composite_ode_interpolant(
        Θ, integrator, caches::T, current, idxs,
        deriv
    ) where {T <: Tuple}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(
            expr.args,
            quote
                if $i == current
                    return ode_interpolant(
                        Θ, integrator.dt, integrator.uprev,
                        integrator.u, integrator.k, caches[$i], idxs,
                        deriv, integrator.differential_vars
                    )
                end
            end
        )
    end
    push!(
        expr.args,
        quote
            throw("Cache $current is not available. There are only $(length(caches)) caches.")
        end
    )
    return expr
end

@inline function ode_interpolant!(val, Θ, integrator::SciMLBase.DEIntegrator, idxs, deriv)
    SciMLBase.addsteps!(integrator)
    return if integrator.cache isa CompositeCache
        ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev, integrator.u,
            integrator.k, integrator.cache.caches[integrator.cache.current],
            idxs, deriv, integrator.differential_vars
        )
    elseif integrator.cache isa DefaultCache
        alg_choice = integrator.cache.current
        if alg_choice == 1
            ode_interpolant!(
                val, Θ, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache.cache1,
                idxs, deriv, integrator.differential_vars
            )
        elseif alg_choice == 2
            ode_interpolant!(
                val, Θ, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache.cache2,
                idxs, deriv, integrator.differential_vars
            )
        elseif alg_choice == 3
            ode_interpolant!(
                val, Θ, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache.cache3,
                idxs, deriv, integrator.differential_vars
            )
        elseif alg_choice == 4
            ode_interpolant!(
                val, Θ, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache.cache4,
                idxs, deriv, integrator.differential_vars
            )
        elseif alg_choice == 5
            ode_interpolant!(
                val, Θ, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache.cache5,
                idxs, deriv, integrator.differential_vars
            )
        elseif alg_choice == 6
            ode_interpolant!(
                val, Θ, integrator.dt, integrator.uprev, integrator.u,
                integrator.k, integrator.cache.cache6,
                idxs, deriv, integrator.differential_vars
            )
        else
            error("DefaultCache invalid alg_choice. File an issue.")
        end
    else
        ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev, integrator.u,
            integrator.k, integrator.cache, idxs, deriv, integrator.differential_vars
        )
    end
end

function default_ode_interpolant!(
        val, Θ, integrator, cache::DefaultCache, alg_choice, idxs, deriv
    )
    if alg_choice == 1
        return ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache1, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 2
        return ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache2, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 3
        return ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache3, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 4
        return ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache4, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 5
        return ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache5, idxs,
            deriv, integrator.differential_vars
        )
    elseif alg_choice == 6
        return ode_interpolant!(
            val, Θ, integrator.dt, integrator.uprev,
            integrator.u, integrator.k, cache.cache6, idxs,
            deriv, integrator.differential_vars
        )
    else
        error("DefaultCache invalid alg_choice. File an issue.")
    end
end

@generated function composite_ode_interpolant!(
        val, Θ, integrator, caches::T, current, idxs,
        deriv
    ) where {T <: Tuple}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(
            expr.args,
            quote
                if $i == current
                    return ode_interpolant!(
                        val, Θ, integrator.dt, integrator.uprev,
                        integrator.u, integrator.k, caches[$i], idxs,
                        deriv
                    )
                end
            end
        )
    end
    push!(
        expr.args,
        quote
            throw("Cache $current is not available. There are only $(length(caches)) caches.")
        end
    )
    return expr
end

@inline function current_interpolant(
        t::Number, integrator::SciMLBase.DEIntegrator, idxs,
        deriv
    )
    Θ = (t - integrator.tprev) / integrator.dt
    return ode_interpolant(Θ, integrator, idxs, deriv)
end

@inline function current_interpolant(t, integrator::SciMLBase.DEIntegrator, idxs, deriv)
    Θ = (t .- integrator.tprev) ./ integrator.dt
    return [ode_interpolant(ϕ, integrator, idxs, deriv) for ϕ in Θ]
end

@inline function current_interpolant!(
        val, t::Number, integrator::SciMLBase.DEIntegrator,
        idxs, deriv
    )
    Θ = (t - integrator.tprev) / integrator.dt
    return ode_interpolant!(val, Θ, integrator, idxs, deriv)
end

@inline function current_interpolant!(
        val, t, integrator::SciMLBase.DEIntegrator, idxs,
        deriv
    )
    Θ = (t .- integrator.tprev) ./ integrator.dt
    return [ode_interpolant!(val, ϕ, integrator, idxs, deriv) for ϕ in Θ]
end

@inline function current_interpolant!(
        val, t::Array, integrator::SciMLBase.DEIntegrator,
        idxs, deriv
    )
    Θ = similar(t)
    @inbounds @simd ivdep for i in eachindex(t)
        Θ[i] = (t[i] - integrator.tprev) / integrator.dt
    end
    return [ode_interpolant!(val, ϕ, integrator, idxs, deriv) for ϕ in Θ]
end

@inline function current_extrapolant(
        t::Number, integrator::SciMLBase.DEIntegrator,
        idxs = nothing, deriv = Val{0}
    )
    Θ = (t - integrator.tprev) / (integrator.t - integrator.tprev)
    return ode_extrapolant(Θ, integrator, idxs, deriv)
end

@inline function current_extrapolant!(
        val, t::Number, integrator::SciMLBase.DEIntegrator,
        idxs = nothing, deriv = Val{0}
    )
    Θ = (t - integrator.tprev) / (integrator.t - integrator.tprev)
    return ode_extrapolant!(val, Θ, integrator, idxs, deriv)
end

@inline function current_extrapolant(
        t::AbstractArray, integrator::SciMLBase.DEIntegrator,
        idxs = nothing, deriv = Val{0}
    )
    Θ = (t .- integrator.tprev) ./ (integrator.t - integrator.tprev)
    return [ode_extrapolant(ϕ, integrator, idxs, deriv) for ϕ in Θ]
end

@inline function current_extrapolant!(
        val, t, integrator::SciMLBase.DEIntegrator,
        idxs = nothing, deriv = Val{0}
    )
    Θ = (t .- integrator.tprev) ./ (integrator.t - integrator.tprev)
    return [ode_extrapolant!(val, ϕ, integrator, idxs, deriv) for ϕ in Θ]
end

@inline function ode_extrapolant!(val, Θ, integrator::SciMLBase.DEIntegrator, idxs, deriv)
    SciMLBase.addsteps!(integrator)
    return if integrator.cache isa CompositeCache
        composite_ode_extrapolant!(
            val, Θ, integrator, integrator.cache.caches,
            integrator.cache.current, idxs, deriv
        )
    elseif integrator.cache isa DefaultCache
        default_ode_extrapolant!(
            val, Θ, integrator, integrator.cache,
            integrator.cache.current, idxs, deriv
        )
    else
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev, integrator.uprev2,
            integrator.uprev, integrator.k, integrator.cache, idxs, deriv, integrator.differential_vars
        )
    end
end

function default_ode_extrapolant!(
        val, Θ, integrator, cache::DefaultCache, alg_choice, idxs, deriv
    )
    return if alg_choice == 1
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache1, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 2
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache2, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 3
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache3, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 4
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache4, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 5
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache5, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 6
        ode_interpolant!(
            val, Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache6, idxs, deriv, integrator.differential_vars
        )
    else
        error("DefaultCache invalid alg_choice. File an issue.")
    end
end

@generated function composite_ode_extrapolant!(
        val, Θ, integrator, caches::T, current, idxs,
        deriv
    ) where {T <: Tuple}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(
            expr.args,
            quote
                if $i == current
                    return ode_interpolant!(
                        val, Θ, integrator.t - integrator.tprev,
                        integrator.uprev2, integrator.uprev,
                        integrator.k, caches[$i], idxs, deriv, integrator.differential_vars
                    )
                end
            end
        )
    end
    push!(
        expr.args,
        quote
            throw("Cache $current is not available. There are only $(length(caches)) caches.")
        end
    )
    return expr
end

@inline function ode_extrapolant(Θ, integrator::SciMLBase.DEIntegrator, idxs, deriv)
    SciMLBase.addsteps!(integrator)
    return if integrator.cache isa CompositeCache
        composite_ode_extrapolant(
            Θ, integrator, integrator.cache.caches,
            integrator.cache.current, idxs, deriv
        )
    elseif integrator.cache isa DefaultCache
        default_ode_extrapolant(
            Θ, integrator, integrator.cache,
            integrator.cache.current, idxs, deriv
        )
    else
        ode_interpolant(
            Θ, integrator.t - integrator.tprev, integrator.uprev2,
            integrator.uprev, integrator.k, integrator.cache, idxs, deriv, integrator.differential_vars
        )
    end
end

function default_ode_extrapolant(
        Θ, integrator, cache::DefaultCache, alg_choice, idxs, deriv
    )
    return if alg_choice == 1
        ode_interpolant(
            Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache1, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 2
        ode_interpolant(
            Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache2, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 3
        ode_interpolant(
            Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache3, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 4
        ode_interpolant(
            Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache4, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 5
        ode_interpolant(
            Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache5, idxs, deriv, integrator.differential_vars
        )
    elseif alg_choice == 6
        ode_interpolant(
            Θ, integrator.t - integrator.tprev,
            integrator.uprev2, integrator.uprev,
            integrator.k, cache.cache6, idxs, deriv, integrator.differential_vars
        )
    else
        error("DefaultCache invalid alg_choice. File an issue.")
    end
end

@generated function composite_ode_extrapolant(
        Θ, integrator, caches::T, current, idxs,
        deriv
    ) where {T <: Tuple}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(
            expr.args,
            quote
                if $i == current
                    return ode_interpolant(
                        Θ, integrator.t - integrator.tprev,
                        integrator.uprev2, integrator.uprev,
                        integrator.k, caches[$i], idxs, deriv, integrator.differential_vars
                    )
                end
            end
        )
    end
    push!(
        expr.args,
        quote
            throw("Cache $current is not available. There are only $(length(caches)) caches.")
        end
    )
    return expr
end

function _evaluate_interpolant(
        f::F, Θ, dt, timeseries, i₋, i₊,
        cache, idxs,
        deriv, ks, ts, p, differential_vars
    ) where {F}
    _ode_addsteps!(
        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
        cache
    ) # update the kcurrent
    return ode_interpolant(
        Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
        cache, idxs, deriv, differential_vars
    )
end
function evaluate_composite_cache(
        f::F, Θ, dt, timeseries, i₋, i₊,
        caches::Tuple{C1, C2, Vararg}, idxs,
        deriv, ks, ts, p, cacheid, differential_vars
    ) where {F, C1, C2}
    if (cacheid -= 1) != 0
        return evaluate_composite_cache(
            f, Θ, dt, timeseries, i₋, i₊, Base.tail(caches),
            idxs,
            deriv, ks, ts, p, cacheid, differential_vars
        )
    end
    return _evaluate_interpolant(
        f, Θ, dt, timeseries, i₋, i₊,
        first(caches), idxs,
        deriv, ks, ts, p, differential_vars
    )
end
function evaluate_composite_cache(
        f::F, Θ, dt, timeseries, i₋, i₊,
        caches::Tuple{C}, idxs,
        deriv, ks, ts, p, _, differential_vars
    ) where {F, C}
    return _evaluate_interpolant(
        f, Θ, dt, timeseries, i₋, i₊,
        only(caches), idxs,
        deriv, ks, ts, p, differential_vars
    )
end

function evaluate_default_cache(
        f::F, Θ, dt, timeseries, i₋, i₊,
        cache::DefaultCache, idxs, deriv, ks, ts, p, cacheid, differential_vars
    ) where {F}
    if cacheid == 1
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache.cache1, idxs, deriv, ks, ts, p, differential_vars
        )
    elseif cacheid == 2
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache.cache2, idxs, deriv, ks, ts, p, differential_vars
        )
    elseif cacheid == 3
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache.cache3, idxs, deriv, ks, ts, p, differential_vars
        )
    elseif cacheid == 4
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache.cache4, idxs, deriv, ks, ts, p, differential_vars
        )
    elseif cacheid == 5
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache.cache5, idxs, deriv, ks, ts, p, differential_vars
        )
    elseif cacheid == 6
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache.cache6, idxs, deriv, ks, ts, p, differential_vars
        )
    end
end

function evaluate_interpolant(
        f::F, Θ, dt, timeseries, i₋, i₊, cache, idxs,
        deriv, ks, ts, id, p, differential_vars
    ) where {F}
    if isdiscretecache(cache)
        return ode_interpolant(
            Θ, dt, timeseries[i₋], timeseries[i₊], 0, cache, idxs,
            deriv, differential_vars
        )
    elseif !id.dense
        return linear_interpolant(Θ, dt, timeseries[i₋], timeseries[i₊], idxs, deriv)
    elseif cache isa CompositeCache
        return evaluate_composite_cache(
            f, Θ, dt, timeseries, i₋, i₊, cache.caches, idxs,
            deriv, ks, ts, p, id.alg_choice[i₊], differential_vars
        )
    elseif cache isa DefaultCache
        return evaluate_default_cache(
            f, Θ, dt, timeseries, i₋, i₊, cache, idxs,
            deriv, ks, ts, p, id.alg_choice[i₊], differential_vars
        )
    else
        return _evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊,
            cache, idxs, deriv, ks, ts, p, differential_vars
        )
    end
end

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation(
        tvals, id::I, idxs, ::Type{deriv}, p,
        continuity::Symbol = :left
    ) where {I, deriv}
    (; ts, timeseries, ks, f, cache, differential_vars) = id
    @inbounds tdir = sign(ts[end] - ts[1])
    idx = sortperm(tvals, rev = tdir < 0)
    # start the search thinking it's ts[1]-ts[2]
    i₋₊ref = Ref((1, 2))
    vals = map(idx) do j
        t = tvals[j]
        (i₋, i₊) = i₋₊ref[]
        if continuity === :left
            # we have i₋ = i₊ = 1 if t = ts[1], i₊ = i₋ + 1 = lastindex(ts) if t > ts[end],
            # and otherwise i₋ and i₊ satisfy ts[i₋] < t ≤ ts[i₊]
            i₊ = min(lastindex(ts), _searchsortedfirst(ts, t, i₊, tdir > 0))
            i₋ = i₊ > 1 ? i₊ - 1 : i₊
        else
            # we have i₋ = i₊ - 1 = 1 if t < ts[1], i₊ = i₋ = lastindex(ts) if t = ts[end],
            # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ t < ts[i₊]
            i₋ = max(1, _searchsortedlast(ts, t, i₋, tdir > 0))
            i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
        end
        id.sensitivitymode && error(SENSITIVITY_INTERP_MESSAGE)
        i₋₊ref[] = (i₋, i₊)
        dt = ts[i₊] - ts[i₋]
        Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t - ts[i₋]) / dt
        evaluate_interpolant(
            f, Θ, dt, timeseries, i₋, i₊, cache, idxs,
            deriv, ks, ts, id, p, differential_vars
        )
    end
    invpermute!(vals, idx)
    return DiffEqArray(vals, tvals)
end

"""
ode_interpolation(tvals,ts,timeseries,ks)

Get the value at tvals where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation!(
        vals, tvals, id::I, idxs, ::Type{deriv}, p,
        continuity::Symbol = :left
    ) where {I, deriv}
    (; ts, timeseries, ks, f, cache, differential_vars) = id
    @inbounds tdir = sign(ts[end] - ts[1])
    idx = sortperm(tvals, rev = tdir < 0)

    # start the search thinking it's in ts[1]-ts[2]
    i₋ = 1
    i₊ = 2
    # if CompositeCache, have an inplace cache for lower allocations
    # (expecting the same algorithms for large portions of ts)
    current_alg = nothing
    cache_i₊ = nothing
    if cache isa CompositeCache
        current_alg = id.alg_choice[i₊]
        cache_i₊ = cache.caches[current_alg]
    elseif cache isa DefaultCache
        current_alg = id.alg_choice[i₊]
    else
        cache_i₊ = cache
    end
    @inbounds for j in idx
        t = tvals[j]

        if continuity === :left
            # we have i₋ = i₊ = 1 if t = ts[1], i₊ = i₋ + 1 = lastindex(ts) if t > ts[end],
            # and otherwise i₋ and i₊ satisfy ts[i₋] < t ≤ ts[i₊]
            i₊ = min(lastindex(ts), _searchsortedfirst(ts, t, i₊, tdir > 0))
            i₋ = i₊ > 1 ? i₊ - 1 : i₊
        else
            # we have i₋ = i₊ - 1 = 1 if t < ts[1], i₊ = i₋ = lastindex(ts) if t = ts[end],
            # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ t < ts[i₊]
            i₋ = max(1, _searchsortedlast(ts, t, i₋, tdir > 0))
            i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
        end
        id.sensitivitymode && error(SENSITIVITY_INTERP_MESSAGE)

        dt = ts[i₊] - ts[i₋]
        Θ = iszero(dt) ? oneunit(t) / oneunit(dt) : (t - ts[i₋]) / dt

        if isdiscretecache(cache)
            if eltype(vals) <: AbstractArray
                ode_interpolant!(
                    vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], 0, cache,
                    idxs, deriv, differential_vars
                )
            else
                vals[j] = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], 0, cache,
                    idxs, deriv, differential_vars
                )
            end
        elseif !id.dense
            if eltype(vals) <: AbstractArray
                linear_interpolant!(
                    vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], idxs,
                    deriv
                )
            else
                vals[j] = linear_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], idxs,
                    deriv
                )
            end
        elseif cache isa DefaultCache
            if current_alg != id.alg_choice[i₊] # switched algorithm
                current_alg = id.alg_choice[i₊]
                if current_alg == 1
                    _ode_addsteps!(
                        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                        cache.cache1
                    ) # update the kcurrent
                    if eltype(vals) <: AbstractArray
                        ode_interpolant!(
                            vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache1, idxs, deriv, differential_vars
                        )
                    else
                        vals[j] = ode_interpolant(
                            Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache1, idxs, deriv, differential_vars
                        )
                    end
                elseif current_alg == 2
                    _ode_addsteps!(
                        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                        cache.cache2
                    ) # update the kcurrent
                    if eltype(vals) <: AbstractArray
                        ode_interpolant!(
                            vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache2, idxs, deriv, differential_vars
                        )
                    else
                        vals[j] = ode_interpolant(
                            Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache2, idxs, deriv, differential_vars
                        )
                    end
                elseif current_alg == 3
                    _ode_addsteps!(
                        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                        cache.cache3
                    ) # update the kcurrent
                    if eltype(vals) <: AbstractArray
                        ode_interpolant!(
                            vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache3, idxs, deriv, differential_vars
                        )
                    else
                        vals[j] = ode_interpolant(
                            Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache3, idxs, deriv, differential_vars
                        )
                    end
                elseif current_alg == 4
                    _ode_addsteps!(
                        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                        cache.cache4
                    ) # update the kcurrent
                    if eltype(vals) <: AbstractArray
                        ode_interpolant!(
                            vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache4, idxs, deriv, differential_vars
                        )
                    else
                        vals[j] = ode_interpolant(
                            Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache4, idxs, deriv, differential_vars
                        )
                    end
                elseif current_alg == 5
                    _ode_addsteps!(
                        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                        cache.cache5
                    ) # update the kcurrent
                    if eltype(vals) <: AbstractArray
                        ode_interpolant!(
                            vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache5, idxs, deriv, differential_vars
                        )
                    else
                        vals[j] = ode_interpolant(
                            Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache5, idxs, deriv, differential_vars
                        )
                    end
                elseif current_alg == 6
                    _ode_addsteps!(
                        ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                        cache.cache6
                    ) # update the kcurrent
                    if eltype(vals) <: AbstractArray
                        ode_interpolant!(
                            vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache6, idxs, deriv, differential_vars
                        )
                    else
                        vals[j] = ode_interpolant(
                            Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                            cache.cache6, idxs, deriv, differential_vars
                        )
                    end
                end
            end
        else
            if cache isa CompositeCache
                if current_alg != id.alg_choice[i₊] # switched algorithm
                    current_alg = id.alg_choice[i₊]
                    @inbounds cache_i₊ = cache.caches[current_alg] # this alloc is costly
                end
            end
            _ode_addsteps!(
                ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                cache_i₊
            ) # update the kcurrent
            if eltype(vals) <: AbstractArray
                ode_interpolant!(
                    vals[j], Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache_i₊, idxs, deriv, differential_vars
                )
            else
                vals[j] = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache_i₊, idxs, deriv, differential_vars
                )
            end
        end
    end

    return vals
end

"""
ode_interpolation(tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation(
        tval::Number, id::I, idxs, ::Type{deriv}, p,
        continuity::Symbol = :left
    ) where {I, deriv}
    (; ts, timeseries, ks, f, cache, differential_vars) = id
    @inbounds tdir = sign(ts[end] - ts[1])

    if continuity === :left
        # we have i₋ = i₊ = 1 if tval = ts[1], i₊ = i₋ + 1 = lastindex(ts) if tval > ts[end],
        # and otherwise i₋ and i₊ satisfy ts[i₋] < tval ≤ ts[i₊]
        i₊ = min(lastindex(ts), _searchsortedfirst(ts, tval, 2, tdir > 0))
        i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
        # we have i₋ = i₊ - 1 = 1 if tval < ts[1], i₊ = i₋ = lastindex(ts) if tval = ts[end],
        # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ tval < ts[i₊]
        i₋ = max(1, _searchsortedlast(ts, tval, 1, tdir > 0))
        i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end
    id.sensitivitymode && error(SENSITIVITY_INTERP_MESSAGE)

    @inbounds begin
        dt = ts[i₊] - ts[i₋]
        Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval - ts[i₋]) / dt

        if isdiscretecache(cache)
            val = ode_interpolant(
                Θ, dt, timeseries[i₋], timeseries[i₊], 0, cache, idxs,
                deriv, differential_vars
            )
        elseif !id.dense
            val = linear_interpolant(Θ, dt, timeseries[i₋], timeseries[i₊], idxs, deriv)
        elseif cache isa CompositeCache
            _ode_addsteps!(
                ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                cache.caches[id.alg_choice[i₊]]
            ) # update the kcurrent
            val = ode_interpolant(
                Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                cache.caches[id.alg_choice[i₊]], idxs, deriv, differential_vars
            )
        elseif cache isa DefaultCache
            alg_choice = id.alg_choice[i₊]
            if alg_choice == 1
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache1
                ) # update the kcurrent
                val = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache1, idxs, deriv, differential_vars
                )
            elseif alg_choice == 2
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache2
                ) # update the kcurrent
                val = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache2, idxs, deriv, differential_vars
                )
            elseif alg_choice == 3
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache3
                ) # update the kcurrent
                val = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache3, idxs, deriv, differential_vars
                )
            elseif alg_choice == 4
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache4
                ) # update the kcurrent
                val = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache4, idxs, deriv, differential_vars
                )
            elseif alg_choice == 5
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache5
                ) # update the kcurrent
                val = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache5, idxs, deriv, differential_vars
                )
            elseif alg_choice == 6
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache6
                ) # update the kcurrent
                val = ode_interpolant(
                    Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache6, idxs, deriv, differential_vars
                )
            else
                error("DefaultCache invalid alg_choice. File an issue.")
            end
        else
            _ode_addsteps!(
                ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                cache
            ) # update the kcurrent
            val = ode_interpolant(
                Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊], cache,
                idxs, deriv, differential_vars
            )
        end
    end

    return val
end

"""
ode_interpolation!(out,tval::Number,ts,timeseries,ks)

Get the value at tval where the solution is known at the
times ts (sorted), with values timeseries and derivatives ks
"""
function ode_interpolation!(
        out, tval::Number, id::I, idxs, ::Type{deriv}, p,
        continuity::Symbol = :left
    ) where {I, deriv}
    (; ts, timeseries, ks, f, cache, differential_vars) = id
    @inbounds tdir = sign(ts[end] - ts[1])

    if continuity === :left
        # we have i₋ = i₊ = 1 if tval = ts[1], i₊ = i₋ + 1 = lastindex(ts) if tval > ts[end],
        # and otherwise i₋ and i₊ satisfy ts[i₋] < tval ≤ ts[i₊]
        i₊ = min(lastindex(ts), _searchsortedfirst(ts, tval, 2, tdir > 0))
        i₋ = i₊ > 1 ? i₊ - 1 : i₊
    else
        # we have i₋ = i₊ - 1 = 1 if tval < ts[1], i₊ = i₋ = lastindex(ts) if tval = ts[end],
        # and otherwise i₋ and i₊ satisfy ts[i₋] ≤ tval < ts[i₊]
        i₋ = max(1, _searchsortedlast(ts, tval, 1, tdir > 0))
        i₊ = i₋ < lastindex(ts) ? i₋ + 1 : i₋
    end
    id.sensitivitymode && error(SENSITIVITY_INTERP_MESSAGE)

    @inbounds begin
        dt = ts[i₊] - ts[i₋]
        Θ = iszero(dt) ? oneunit(tval) / oneunit(dt) : (tval - ts[i₋]) / dt

        if isdiscretecache(cache)
            ode_interpolant!(
                out, Θ, dt, timeseries[i₋], timeseries[i₊], 0, cache, idxs,
                deriv, differential_vars
            )
        elseif !id.dense
            linear_interpolant!(out, Θ, dt, timeseries[i₋], timeseries[i₊], idxs, deriv)
        elseif cache isa CompositeCache
            _ode_addsteps!(
                ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                cache.caches[id.alg_choice[i₊]]
            ) # update the kcurrent
            ode_interpolant!(
                out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                cache.caches[id.alg_choice[i₊]], idxs, deriv, differential_vars
            )
        elseif cache isa DefaultCache
            alg_choice = id.alg_choice[i₊]
            if alg_choice == 1
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache1
                ) # update the kcurrent
                ode_interpolant!(
                    out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache1, idxs, deriv, differential_vars
                )
            elseif alg_choice == 2
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache2
                ) # update the kcurrent
                ode_interpolant!(
                    out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.caches[2], idxs, deriv, differential_vars
                )
            elseif alg_choice == 3
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache3
                ) # update the kcurrent
                ode_interpolant!(
                    out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache3, idxs, deriv, differential_vars
                )
            elseif alg_choice == 4
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache4
                ) # update the kcurrent
                ode_interpolant!(
                    out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache5, idxs, deriv, differential_vars
                )
            elseif alg_choice == 5
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache5
                ) # update the kcurrent
                ode_interpolant!(
                    out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache5, idxs, deriv, differential_vars
                )
            elseif alg_choice == 6
                _ode_addsteps!(
                    ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                    cache.cache6
                ) # update the kcurrent
                ode_interpolant!(
                    out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊],
                    cache.cache6, idxs, deriv, differential_vars
                )
            else
                error("DefaultCache invalid alg_choice. File an issue.")
            end
        else
            _ode_addsteps!(
                ks[i₊], ts[i₋], timeseries[i₋], timeseries[i₊], dt, f, p,
                cache
            ) # update the kcurrent
            ode_interpolant!(
                out, Θ, dt, timeseries[i₋], timeseries[i₊], ks[i₊], cache,
                idxs, deriv, differential_vars
            )
        end
    end

    return out
end

"""
By default, Hermite interpolant so update the derivative at the two ends
"""
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache, always_calc_begin = false,
        allow_calc_end = true, force_calc_end = false
    )
    if length(k) < 2 || always_calc_begin
        if cache isa OrdinaryDiffEqMutableCache
            rtmp = similar(u, eltype(eltype(k)))
            f(rtmp, uprev, p, t)
            copyat_or_push!(k, 1, rtmp)
            f(rtmp, u, p, t + dt)
            copyat_or_push!(k, 2, rtmp)
        else
            copyat_or_push!(k, 1, f(uprev, p, t))
            copyat_or_push!(k, 2, f(u, p, t + dt))
        end
    end
    return nothing
end

"""
ode_interpolant and ode_interpolant! dispatch
"""
function ode_interpolant(
        Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{TI}}, differential_vars
    ) where {TI}
    return _ode_interpolant(Θ, dt, y₀, y₁, k, cache, idxs, T, differential_vars)
end

function ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::OrdinaryDiffEqMutableCache, idxs,
        T::Type{Val{TI}}, differential_vars
    ) where {TI}
    return if idxs isa Number || y₀ isa Union{Number, SArray}
        # typeof(y₀) can be these if saveidxs gives a single value
        _ode_interpolant(Θ, dt, y₀, y₁, k, cache, idxs, T, differential_vars)
    elseif idxs isa Nothing
        if y₁ isa Array{<:Number}
            out = similar(y₁, eltype(first(y₁) * oneunit(Θ)))
            copyto!(out, y₁)
        else
            out = oneunit(Θ) .* y₁
        end
        _ode_interpolant!(out, Θ, dt, y₀, y₁, k, cache, idxs, T, differential_vars)
    else
        if y₁ isa Array{<:Number}
            out = similar(y₁, eltype(first(y₁) * oneunit(Θ)), axes(idxs))
            for i in eachindex(idxs)
                out[i] = y₁[idxs[i]]
            end
        else
            out = oneunit(Θ) .* y₁[idxs]
        end
        _ode_interpolant!(out, Θ, dt, y₀, y₁, k, cache, idxs, T, differential_vars)
    end
end

function ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{TI}}, differential_vars
    ) where {TI}
    return _ode_interpolant!(out, Θ, dt, y₀, y₁, k, cache, idxs, T, differential_vars)
end

##################### Hermite Interpolants

function interpolation_differential_vars(differential_vars, y₀, idxs)
    if isnothing(differential_vars)
        if y₀ isa Number
            return true
        elseif idxs === nothing
            return Trues(size(y₀))
        elseif idxs isa Number
            return true
        else
            return Trues(size(idxs))
        end
    elseif differential_vars isa DifferentialVarsUndefined #for non diagonal mass matrices, use linear interpolation.
        if y₀ isa Number
            return false
        elseif idxs === nothing
            return Falses(size(y₀))
        elseif idxs isa Number
            return false
        else
            return Falses(size(idxs))
        end
    elseif idxs isa Number
        return return differential_vars[idxs]
    elseif idxs === nothing
        return differential_vars
    else
        return @view differential_vars[idxs]
    end
end

# If no dispatch found, assume Hermite (or linear when k is empty, e.g. SDE)
function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{TI}}, differential_vars
    ) where {TI}
    TI > 3 && throw(DerivativeOrderNotPossibleError())

    # Linear fallback when no dense output vectors (e.g. SDE)
    if isempty(k)
        return linear_interpolant(Θ, dt, y₀, y₁, idxs, T)
    end

    differential_vars = interpolation_differential_vars(differential_vars, y₀, idxs)
    return hermite_interpolant(
        Θ, dt, y₀, y₁, k, Val{cache isa OrdinaryDiffEqMutableCache},
        idxs, T, differential_vars
    )
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{TI}}, differential_vars
    ) where {TI}
    TI > 3 && throw(DerivativeOrderNotPossibleError())

    if isempty(k)
        return linear_interpolant!(out, Θ, dt, y₀, y₁, idxs, T)
    end

    differential_vars = interpolation_differential_vars(differential_vars, y₀, idxs)
    return hermite_interpolant!(out, Θ, dt, y₀, y₁, k, idxs, T, differential_vars)
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{false}}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars
    )
    #@.. broadcast=false (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
    if all(differential_vars)
        @inbounds (1 - Θ) * y₀ + Θ * y₁ +
            (
            Θ * (Θ - 1) *
                ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
        )
    else
        @inbounds (1 - Θ) * y₀ + Θ * y₁ +
            differential_vars .* (
            Θ * (Θ - 1) *
                ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
        )
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{true}}, idxs::Nothing,
        T::Type{Val{0}}, differential_vars
    )
    #@.. broadcast=false (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*k[1] + Θ*dt*k[2])
    if all(differential_vars)
        @inbounds @.. broadcast = false (1 - Θ) * y₀ + Θ * y₁ +
            Θ * (Θ - 1) *
            ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
    else
        @inbounds @.. broadcast = false (1 - Θ) * y₀ + Θ * y₁ +
            differential_vars * Θ * (Θ - 1) *
            ((1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] + Θ * dt * k[2])
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀::Array, y₁, k, ::Type{Val{true}},
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    out = similar(y₀)
    @inbounds @simd ivdep for i in eachindex(y₀)
        out[i] = (1 - Θ) * y₀[i] + Θ * y₁[i] +
            differential_vars[i] * Θ * (Θ - 1) *
            ((1 - 2Θ) * (y₁[i] - y₀[i]) + (Θ - 1) * dt * k[1][i] + Θ * dt * k[2][i])
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{0}}, differential_vars
    )
    # return @.. broadcast=false (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+(Θ-1)*dt*k[1][idxs] + Θ*dt*k[2][idxs])
    if all(differential_vars)
        return (1 - Θ) * y₀[idxs] + Θ * y₁[idxs] +
            (
            Θ * (Θ - 1) *
                (
                (1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] +
                    Θ * dt * k[2][idxs]
            )
        )
    else
        return (1 - Θ) * y₀[idxs] + Θ * y₁[idxs] +
            differential_vars .* (
            Θ * (Θ - 1) *
                (
                (1 - 2Θ) * (y₁[idxs] - y₀[idxs]) + (Θ - 1) * dt * k[1][idxs] +
                    Θ * dt * k[2][idxs]
            )
        )
    end
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false out = (1 - Θ) * y₀ + Θ * y₁ +
            Θ * (Θ - 1) *
            (
            (1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] +
                Θ * dt * k[2]
        )
    else
        @inbounds @.. broadcast = false out = (1 - Θ) * y₀ + Θ * y₁ +
            differential_vars * Θ * (Θ - 1) *
            (
            (1 - 2Θ) * (y₁ - y₀) + (Θ - 1) * dt * k[1] +
                Θ * dt * k[2]
        )
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs::Nothing,
        T::Type{Val{0}}, differential_vars
    )
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = (1 - Θ) * y₀[i] + Θ * y₁[i] +
            differential_vars[i] * Θ * (Θ - 1) *
            ((1 - 2Θ) * (y₁[i] - y₀[i]) + (Θ - 1) * dt * k[1][i] + Θ * dt * k[2][i])
    end
    out
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{0}}, differential_vars
    )
    if all(differential_vars)
        @views @.. broadcast = false out = (1 - Θ) * y₀[idxs] + Θ * y₁[idxs] +
            Θ * (Θ - 1) *
            (
            (1 - 2Θ) * (y₁[idxs] - y₀[idxs]) +
                (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs]
        )
    else
        @views @.. broadcast = false out = (1 - Θ) * y₀[idxs] + Θ * y₁[idxs] +
            differential_vars * Θ * (Θ - 1) *
            (
            (1 - 2Θ) * (y₁[idxs] - y₀[idxs]) +
                (Θ - 1) * dt * k[1][idxs] + Θ * dt * k[2][idxs]
        )
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{0}}, differential_vars
    )
    @inbounds for (j, i) in enumerate(idxs)
        out[j] = (1 - Θ) * y₀[i] + Θ * y₁[i] +
            differential_vars[j] * Θ * (Θ - 1) *
            ((1 - 2Θ) * (y₁[i] - y₀[i]) + (Θ - 1) * dt * k[1][i] + Θ * dt * k[2][i])
    end
    out
end

"""
Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{false}}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars
    )
    #@.. broadcast=false k[1] + Θ*(-4*dt*k[1] - 2*dt*k[2] - 6*y₀ + Θ*(3*dt*k[1] + 3*dt*k[2] + 6*y₀ - 6*y₁) + 6*y₁)/dt
    if all(differential_vars)
        @inbounds (
            k[1] +
                Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ * (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) + 6 * y₁
            ) / dt
        )
    else
        @inbounds (.!differential_vars) .* (y₁ - y₀) / dt +
            differential_vars .* (
            k[1] +
                Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ * (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) + 6 * y₁
            ) / dt
        )
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{true}}, idxs::Nothing,
        T::Type{Val{1}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false !differential_vars *
            (
            (y₁ - y₀) /
                dt
        ) + (
            k[1] +
                Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ *
                    (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                    6 * y₁
            ) / dt
        )
    else
        @inbounds @.. broadcast = false !differential_vars *
            (
            (y₁ - y₀) /
                dt
        ) + differential_vars * (
            k[1] +
                Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ *
                    (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                    6 * y₁
            ) / dt
        )
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{1}}, differential_vars
    )
    if all(differential_vars)
        (
            k[1][idxs] +
                Θ * (
                -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] +
                    Θ * (3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] + 6 * y₀[idxs] - 6 * y₁[idxs]) +
                    6 * y₁[idxs]
            ) / dt
        )
    else
        (.!differential_vars) .* ((y₁[idxs] - y₀[idxs]) / dt) +
            differential_vars .* (
            k[1][idxs] +
                Θ * (
                -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] +
                    Θ * (3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] + 6 * y₀[idxs] - 6 * y₁[idxs]) +
                    6 * y₁[idxs]
            ) / dt
        )
    end
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false out = (
            k[1] +
                Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ *
                    (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                    6 * y₁
            ) / dt
        )
    else
        @inbounds @.. broadcast = false out = !differential_vars * ((y₁ - y₀) / dt) +
            differential_vars * (
            k[1] +
                Θ * (
                -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                    Θ *
                    (3 * dt * k[1] + 3 * dt * k[2] + 6 * y₀ - 6 * y₁) +
                    6 * y₁
            ) / dt
        )
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs::Nothing,
        T::Type{Val{1}}, differential_vars
    )
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = !differential_vars[i] * ((y₁[i] - y₀[i]) / dt) +
            differential_vars[i] * (
            k[1][i] +
                Θ * (
                -4 * dt * k[1][i] - 2 * dt * k[2][i] - 6 * y₀[i] +
                    Θ * (3 * dt * k[1][i] + 3 * dt * k[2][i] + 6 * y₀[i] - 6 * y₁[i]) +
                    6 * y₁[i]
            ) / dt
        )
    end
    out
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{1}}, differential_vars
    )
    if all(differential_vars)
        @views @.. broadcast = false out = (
            k[1][idxs] +
                Θ * (
                -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] -
                    6 * y₀[idxs] +
                    Θ * (
                    3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] +
                        6 * y₀[idxs] - 6 * y₁[idxs]
                ) + 6 * y₁[idxs]
            ) / dt
        )
    else
        @views @.. broadcast = false out = !differential_vars * ((y₁ - y₀) / dt) +
            differential_vars * (
            k[1][idxs] +
                Θ * (
                -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] -
                    6 * y₀[idxs] +
                    Θ * (
                    3 * dt * k[1][idxs] + 3 * dt * k[2][idxs] +
                        6 * y₀[idxs] - 6 * y₁[idxs]
                ) + 6 * y₁[idxs]
            ) / dt
        )
    end
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{1}}, differential_vars
    )
    @inbounds for (j, i) in enumerate(idxs)
        out[j] = !differential_vars[j] * ((y₁[i] - y₀[i]) / dt) +
            differential_vars[j] * (
            k[1][i] +
                Θ * (
                -4 * dt * k[1][i] - 2 * dt * k[2][i] - 6 * y₀[i] +
                    Θ * (3 * dt * k[1][i] + 3 * dt * k[2][i] + 6 * y₀[i] - 6 * y₁[i]) +
                    6 * y₁[i]
            ) / dt
        )
    end
    out
end

"""
Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{false}}, idxs::Nothing,
        T::Type{Val{2}}, differential_vars
    )
    if all(differential_vars)
        @inbounds (
            -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                Θ * (6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ - 12 * y₁) + 6 * y₁
        ) /
            (dt * dt)
    else
        @inbounds differential_vars .* (
            -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                Θ * (6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ - 12 * y₁) + 6 * y₁
        ) /
            (dt * dt)
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{true}}, idxs::Nothing,
        T::Type{Val{2}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false (
            -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                Θ *
                (6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ - 12 * y₁) +
                6 * y₁
        ) / (dt * dt)
    else
        @inbounds @.. broadcast = false differential_vars *
            (
            -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                Θ *
                (6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ - 12 * y₁) +
                6 * y₁
        ) / (dt * dt)
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{2}}, differential_vars
    )
    if all(differential_vars)
        @views out = (
            -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] +
                Θ * (
                6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] + 12 * y₀[idxs] -
                    12 * y₁[idxs]
            ) + 6 * y₁[idxs]
        ) / (dt * dt)
    else
        @views out = differential_vars .* (-4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] - 6 * y₀[idxs] + Θ * (6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] + 12 * y₀[idxs] - 12 * y₁[idxs]) + 6 * y₁[idxs]) / (dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs::Nothing, T::Type{Val{2}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false out = (
            -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                Θ *
                (
                6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ -
                    12 * y₁
            ) +
                6 * y₁
        ) / (dt * dt)
    else
        @inbounds @.. broadcast = false out = differential_vars *
            (
            -4 * dt * k[1] - 2 * dt * k[2] - 6 * y₀ +
                Θ *
                (
                6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ -
                    12 * y₁
            ) +
                6 * y₁
        ) / (dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs::Nothing,
        T::Type{Val{2}}, differential_vars
    )
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = differential_vars[i] *
            (
            -4 * dt * k[1][i] - 2 * dt * k[2][i] - 6 * y₀[i] +
                Θ * (6 * dt * k[1][i] + 6 * dt * k[2][i] + 12 * y₀[i] - 12 * y₁[i]) +
                6 * y₁[i]
        ) / (dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{2}}, differential_vars
    )
    if all(differential_vars)
        @views @.. broadcast = false out = (
            -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] -
                6 * y₀[idxs] +
                Θ * (
                6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] +
                    12 * y₀[idxs] - 12 * y₁[idxs]
            ) + 6 * y₁[idxs]
        ) /
            (dt * dt)
    else
        @views @.. broadcast = false out = differential_vars *
            (
            -4 * dt * k[1][idxs] - 2 * dt * k[2][idxs] -
                6 * y₀[idxs] +
                Θ * (
                6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] +
                    12 * y₀[idxs] - 12 * y₁[idxs]
            ) + 6 * y₁[idxs]
        ) /
            (dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{2}}, differential_vars
    )
    @inbounds for (j, i) in enumerate(idxs)
        out[j] = differential_vars[j] *
            (
            -4 * dt * k[1][i] - 2 * dt * k[2][i] - 6 * y₀[i] +
                Θ * (6 * dt * k[1][i] + 6 * dt * k[2][i] + 12 * y₀[i] - 12 * y₁[i]) +
                6 * y₁[i]
        ) / (dt * dt)
    end
    out
end

"""
Herimte Interpolation, chosen if no other dispatch for ode_interpolant
"""
@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{false}}, idxs::Nothing,
        T::Type{Val{3}}, differential_vars
    )
    #@.. broadcast=false (6*dt*k[1] + 6*dt*k[2] + 12*y₀ - 12*y₁)/(dt*dt*dt)
    if all(differential_vars)
        @inbounds (6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ - 12 * y₁) / (dt * dt * dt)
    else
        @inbounds differential_vars .* (6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ - 12 * y₁) /
            (dt * dt * dt)
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, ::Type{Val{true}}, idxs::Nothing,
        T::Type{Val{3}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false (
            6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ -
                12 * y₁
        ) / (
            dt *
                dt *
                dt
        )
    else
        @inbounds @.. broadcast = false differential_vars *
            (
            6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ -
                12 * y₁
        ) / (
            dt *
                dt *
                dt
        )
    end
end

@muladd function hermite_interpolant(
        Θ, dt, y₀, y₁, k, cache, idxs, T::Type{Val{3}}, differential_vars
    )
    if all(differential_vars)
        @views out = (
            6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] + 12 * y₀[idxs] -
                12 * y₁[idxs]
        ) /
            (dt * dt * dt)
    else
        @views out = differential_vars .* (6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] + 12 * y₀[idxs] - 12 * y₁[idxs]) / (dt * dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs::Nothing, T::Type{Val{3}}, differential_vars
    )
    if all(differential_vars)
        @inbounds @.. broadcast = false out = (
            6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ -
                12 * y₁
        ) /
            (dt * dt * dt)
    else
        @inbounds @.. broadcast = false out = differential_vars *
            (
            6 * dt * k[1] + 6 * dt * k[2] + 12 * y₀ -
                12 * y₁
        ) /
            (dt * dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs::Nothing,
        T::Type{Val{3}}, differential_vars
    )
    @inbounds @simd ivdep for i in eachindex(out)
        out[i] = differential_vars[i] *
            (6 * dt * k[1][i] + 6 * dt * k[2][i] + 12 * y₀[i] - 12 * y₁[i]) /
            (dt * dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{3}}, differential_vars
    )
    if all(differential_vars)
        @views @.. broadcast = false out = (
            6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] +
                12 * y₀[idxs] - 12 * y₁[idxs]
        ) / (dt * dt * dt)
    else
        @views @.. broadcast = false out = differential_vars *
            (
            6 * dt * k[1][idxs] + 6 * dt * k[2][idxs] +
                12 * y₀[idxs] - 12 * y₁[idxs]
        ) / (dt * dt * dt)
    end
    out
end

@muladd function hermite_interpolant!(
        out::Array, Θ, dt, y₀, y₁, k, idxs, T::Type{Val{3}}, differential_vars
    )
    @inbounds for (j, i) in enumerate(idxs)
        out[j] = differential_vars[j] *
            (6 * dt * k[1][i] + 6 * dt * k[2][i] + 12 * y₀[i] - 12 * y₁[i]) /
            (dt * dt * dt)
    end
    out
end

######################## Linear Interpolants

@muladd @inline function linear_interpolant(Θ, dt, y₀, y₁, idxs::Nothing, T::Type{Val{0}})
    Θm1 = (1 - Θ)
    @.. broadcast = false Θm1 * y₀ + Θ * y₁
end

@muladd @inline function linear_interpolant(Θ, dt, y₀, y₁, idxs, T::Type{Val{0}})
    Θm1 = (1 - Θ)
    @.. broadcast = false Θm1 * y₀[idxs] + Θ * y₁[idxs]
end

@muladd @inline function linear_interpolant!(
        out, Θ, dt, y₀, y₁, idxs::Nothing,
        T::Type{Val{0}}
    )
    Θm1 = (1 - Θ)
    @.. broadcast = false out = Θm1 * y₀ + Θ * y₁
    out
end

@muladd @inline function linear_interpolant!(out, Θ, dt, y₀, y₁, idxs, T::Type{Val{0}})
    Θm1 = (1 - Θ)
    @views @.. broadcast = false out = Θm1 * y₀[idxs] + Θ * y₁[idxs]
    out
end

"""
Linear Interpolation
"""
@inline function linear_interpolant(Θ, dt, y₀, y₁, idxs::Nothing, T::Type{Val{1}})
    return (y₁ - y₀) / dt
end

@inline function linear_interpolant(Θ, dt, y₀, y₁, idxs, T::Type{Val{1}})
    return @.. broadcast = false (y₁[idxs] - y₀[idxs]) / dt
end

@inline function linear_interpolant!(out, Θ, dt, y₀, y₁, idxs::Nothing, T::Type{Val{1}})
    @.. broadcast = false out = (y₁ - y₀) / dt
    return out
end

@inline function linear_interpolant!(out, Θ, dt, y₀, y₁, idxs, T::Type{Val{1}})
    @views @.. broadcast = false out = (y₁[idxs] - y₀[idxs]) / dt
    return out
end
