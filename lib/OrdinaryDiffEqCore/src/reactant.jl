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

function _traced_finalize_solution(integrator, retcode)
    return ConstructionBase.setproperties(
        integrator.sol,
        (;
            u = [integrator.u], t = [integrator.t], k = nothing, prob = nothing,
            interp = nothing, dense = false, stats = nothing, retcode,
        )
    )
end
