_maybe_traced(x) = ReactantCore.within_compile() ? ReactantCore.promote_to_traced(x) : x

# Function wrappers hide types from Reactant and provide no compile-time reuse inside `@jit`.
SciMLBase.specialization(::ODEFunction{iip, SciMLBase.AutoSpecialize}) where {iip} =
    ReactantCore.within_compile() ? SciMLBase.FullSpecialize : SciMLBase.AutoSpecialize

# Reactant requires every loop-carried path to own its traced value at the loop boundary.
function _dealias_traced!(x)
    ReactantCore.within_compile() || return x
    x isa Union{Type, Function, Module, AbstractString, Symbol, SciMLBase.AbstractSciMLFunction} && return x
    if x isa DenseArray
        if eltype(x) <: Number
            return x .* one(eltype(x))
        end
        y = copy(x)
        for i in eachindex(x)
            y[i] = _dealias_traced!(x[i])
        end
        return y
    end
    if x isa Number
        return x * one(x)
    end
    x isa Union{Tuple, NamedTuple} && return map(_dealias_traced!, x)
    T = typeof(x)
    isbitstype(T) && return x
    if ismutable(x)
        for name in fieldnames(T)
            (isdefined(x, name) && !isconst(T, name)) || continue
            setfield!(x, name, _dealias_traced!(getfield(x, name)))
        end
        return x
    end
    names = fieldnames(T)
    isempty(names) && return x
    values = map(name -> _dealias_traced!(getfield(x, name)), names)
    return ConstructionBase.setproperties(x, NamedTuple{names}(values))
end

_traced_select(cond::Bool, a, b) = ifelse(cond, a, b)
_traced_select(cond, a, b) = a isa AbstractArray ? ifelse.(cond, a, b) : ifelse(cond, a, b)

function _traced_update_fsal!(integrator, accepted = true)
    if isinplace(integrator.sol.prob)
        fsalfirst, fsallast = get_fsalfirstlast(integrator.cache, integrator.u)
        if !isnothing(fsalfirst) && !isnothing(fsallast)
            recursivecopy!(fsalfirst, _traced_select(accepted, fsallast, fsalfirst))
        end
    else
        integrator.fsalfirst = _traced_select(
            accepted, integrator.fsallast, integrator.fsalfirst
        )
    end
    return nothing
end

function _traced_finalize_solution(integrator, retcode)
    return ConstructionBase.setproperties(
        integrator.sol,
        (;
            u = [integrator.u], t = [integrator.t], k = nothing, prob = nothing,
            interp = nothing, dense = false, stats = nothing, retcode,
        )
    )
end

function _traced_fixed_step_solve!(integrator)
    tstop = integrator.tdir * first_tstop(integrator)
    nsteps = _maybe_traced(0)
    _dealias_traced!(integrator)
    ReactantCore.@trace track_numbers = false while (
            (integrator.tdir * integrator.t < tstop) &
                (nsteps < integrator.opts.maxiters)
        )
        integrator.dt = integrator.tdir * min(
            abs(integrator.dtcache), tstop - integrator.tdir * integrator.t
        )
        perform_step!(integrator, integrator.cache)
        integrator.tprev = integrator.t
        integrator.t += integrator.dt
        update_uprev!(integrator)
        _traced_update_fsal!(integrator)
        nsteps += one(nsteps)
        _dealias_traced!(integrator)
    end
    retcode = ifelse(
        integrator.tdir * integrator.t >= tstop,
        ReturnCode.Success,
        ReturnCode.MaxIters
    )
    return _traced_finalize_solution(integrator, retcode)
end

function _traced_adaptive_solve!(integrator, cache::PIControllerCache)
    tstop = integrator.tdir * first_tstop(integrator)
    nsteps = _maybe_traced(0)
    naccept = _maybe_traced(0)
    _dealias_traced!(integrator)
    ReactantCore.@trace track_numbers = false while (
            (integrator.tdir * integrator.t < tstop) &
                (nsteps < integrator.opts.maxiters)
        )
        integrator.dt = integrator.tdir * min(
            abs(integrator.dt), tstop - integrator.tdir * integrator.t
        )
        old_fsalfirst = integrator.fsalfirst
        perform_step!(integrator, integrator.cache)

        controller = cache.controller
        (; qmin, qmax, qmax_first_step, gamma, qsteady_min, qsteady_max) =
            controller.basic
        qmax = ifelse(iszero(naccept), qmax_first_step, qmax)
        EEst = SciMLBase.value(get_EEst(integrator))
        q11 = fastpower(EEst, controller.beta1)
        q = q11 / fastpower(cache.errold, controller.beta2)
        q = clamp(q / gamma, inv(qmax), inv(qmin))
        q = ifelse(iszero(EEst), inv(qmax), q)
        accepted = EEst <= one(EEst)

        accepted_q = ifelse((qsteady_min <= q) & (q <= qsteady_max), one(q), q)
        accepted_dt = integrator.dt / accepted_q
        rejected_dt = integrator.dt / min(inv(qmin), q11 / gamma)
        next_dt = ifelse(accepted, accepted_dt, rejected_dt)
        next_dt = integrator.tdir * min(abs(integrator.opts.dtmax), abs(next_dt))
        next_dt = integrator.tdir * max(abs(next_dt), abs(integrator.opts.dtmin))

        cache.q11 = q11
        cache.errold = ifelse(
            accepted, max(EEst, controller.qoldinit), cache.errold
        )
        accepted_u = _traced_select(accepted, integrator.u, integrator.uprev)
        integrator.tprev = ifelse(accepted, integrator.t, integrator.tprev)
        integrator.t = ifelse(accepted, integrator.t + integrator.dt, integrator.t)
        integrator.u = accepted_u
        integrator.uprev = accepted_u
        integrator.fsalfirst = old_fsalfirst
        _traced_update_fsal!(integrator, accepted)
        integrator.dt = next_dt
        nsteps += one(nsteps)
        naccept += ifelse(accepted, one(naccept), zero(naccept))
        _dealias_traced!(integrator)
    end
    retcode = ifelse(
        integrator.tdir * integrator.t >= tstop,
        ReturnCode.Success,
        ReturnCode.MaxIters
    )
    return _traced_finalize_solution(integrator, retcode)
end

function _traced_adaptive_solve!(integrator, cache::AbstractControllerCache)
    throw(
        ArgumentError(
            "$(nameof(typeof(cache))) is not supported inside Reactant compilation; use a PIController"
        )
    )
end

function _traced_ode_initdt_oop(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm, prob, order, integrator
    )
    T = eltype(t)
    oneunit_t = oneunit(t)
    dtmin = max(integrator.opts.dtmin, oneunit_t * eps(T))
    smalldt = max(dtmin, oneunit_t * T(1.0e-6))
    sk = @.. broadcast = false abstol + internalnorm(u0, t) * reltol
    d₀ = internalnorm(u0 ./ sk, t)
    f₀ = prob.f(u0, prob.p, t)
    d₁ = internalnorm(f₀ ./ sk .* oneunit_t, t)
    dt₀ = ifelse(
        (d₀ < T(1.0e-5)) | (d₁ < T(1.0e-5)),
        smalldt,
        oneunit_t * SciMLBase.value(d₀ / d₁) / T(100)
    )
    dt₀ = min(dt₀, tdir * dtmax)
    u₁ = @.. broadcast = false u0 + tdir * dt₀ * f₀
    f₁ = prob.f(u₁, prob.p, t + tdir * dt₀)
    d₂ = internalnorm((f₁ .- f₀) ./ sk .* oneunit_t, t) / dt₀ * oneunit_t
    max_d₁d₂ = max(d₁, d₂)
    dt₁ = ifelse(
        max_d₁d₂ <= T(1.0e-15),
        max(smalldt, dt₀ * T(1.0e-3)),
        oneunit_t * T(10)^(-(2 + log10(max_d₁d₂)) / order)
    )
    return tdir * max(dtmin, min(100dt₀, dt₁, tdir * dtmax))
end

function _traced_ode_initdt_iip(
        u0, t, tdir, dtmax, abstol, reltol, internalnorm, prob, order, integrator
    )
    T = eltype(t)
    oneunit_t = oneunit(t)
    dtmin = max(integrator.opts.dtmin, oneunit_t * eps(T))
    smalldt = max(dtmin, oneunit_t * T(1.0e-6))
    sk = @.. broadcast = false abstol + internalnorm(u0, t) * reltol
    d₀ = internalnorm(u0 ./ sk, t)
    f₀ = zero(u0)
    prob.f(f₀, u0, prob.p, t)
    d₁ = internalnorm(f₀ ./ sk .* oneunit_t, t)
    dt₀ = ifelse(
        (d₀ < T(1.0e-5)) | (d₁ < T(1.0e-5)),
        smalldt,
        oneunit_t * SciMLBase.value(d₀ / d₁) / T(100)
    )
    dt₀ = min(dt₀, tdir * dtmax)
    u₁ = @.. broadcast = false u0 + tdir * dt₀ * f₀
    f₁ = zero(f₀)
    prob.f(f₁, u₁, prob.p, t + tdir * dt₀)
    d₂ = internalnorm((f₁ .- f₀) ./ sk .* oneunit_t, t) / dt₀ * oneunit_t
    max_d₁d₂ = max(d₁, d₂)
    dt₁ = ifelse(
        max_d₁d₂ <= T(1.0e-15),
        max(smalldt, dt₀ * T(1.0e-3)),
        oneunit_t * T(10)^(-(2 + log10(max_d₁d₂)) / order)
    )
    return tdir * max(dtmin, min(100dt₀, dt₁, tdir * dtmax))
end
