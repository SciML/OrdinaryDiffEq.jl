function _check_traced_time_grid(prob, dt, dtmax, adaptive, stage_limiter!)
    tstart, tend = prob.tspan
    tmag = max(abs(tstart), abs(tend))
    tmag isa Union{Float16, Float32, Float64} && issubnormal(eps(tmag)) || return
    span = abs(tend - tstart)
    abs(dtmax) < span && throw(ArgumentError("timespans with floatmin-scale spacing and a bounded dtmax are not supported under Reactant compilation because the backend flushes subnormals; rescale tspan/u0"))
    adaptive && return
    dt isa Real && abs(dt) + tstop_tolerance(tmag) >= span && return
    stage_limiter! === trivial_limiter! && _rhs_ignores_state_and_time(prob) && return
    throw(ArgumentError("Reactant cannot preserve the requested fixed time grid at subnormal time spacing; rescale tspan"))
end

function _rhs_ignores_state_and_time(prob)
    prob.f isa SciMLBase.ODEFunction || return false
    rhs = SciMLBase.unwrapped_f(prob.f.f)
    candidates = methods(rhs)
    length(candidates) == 1 || return false
    method = only(candidates)
    iip = isinplace(prob)
    (method.isva || method.nargs != 4 + iip || isdefined(method, :generator)) && return false
    # In lowered Julia IR, slot 1 is the callable and flag 0x08 marks an argument as used.
    flags = only(code_lowered(rhs)).slotflags
    return iszero((flags[2 + iip] | flags[4 + iip]) & 0x08)
end

_maybe_traced(x) = ReactantCore.within_compile() ? ReactantCore.promote_to_traced(x) : x

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

function _traced_finalize_solution(integrator, retcode)
    return ConstructionBase.setproperties(
        integrator.sol,
        (;
            u = [integrator.u], t = [integrator.t], k = nothing, prob = nothing,
            interp = nothing, dense = false, stats = nothing, retcode,
        )
    )
end
