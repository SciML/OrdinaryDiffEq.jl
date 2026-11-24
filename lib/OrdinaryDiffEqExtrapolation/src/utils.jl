@inline _unwrap_val(::Val{B}) where {B} = B
@inline _unwrap_val(value) = value

@inline _vec(value::Union{Number, AbstractVector, SciMLOperators.AbstractSciMLScalarOperator}) =
    value
@inline _vec(value) = vec(value)

@inline _reshape(value::Union{Number, SciMLOperators.AbstractSciMLScalarOperator}, _) = value
@inline _reshape(value, dimensions) = reshape(value, dimensions)

@inline _has_Wfact(f) = hasproperty(f, :Wfact) && getproperty(f, :Wfact) !== nothing

@inline _float_controller_type(::Type{T}) where {T <: Integer} = float(T)
@inline _float_controller_type(::Type{T}) where {T} = T
@inline _controller_scalar_type(
    ::OrdinaryDiffEqCore.CommonControllerOptions{T}
) where {T} = _float_controller_type(T)
@inline _controller_scalar_type(::NamedTuple{(), Tuple{}}) = Float64
@inline _controller_scalar_type(overrides::NamedTuple) =
    _float_controller_type(promote_type(map(typeof, values(overrides))...))

# Splits `range` the way `@threads :static` does; chunk `c` owns scratch slot `c`.
@inline _nchunks(range) = min(length(range), Threads.nthreads())

@inline function _chunk(range, c, nchunks)
    len, rem = divrem(length(range), nchunks)
    first = (c - 1) * len + min(c - 1, rem) + 1
    return range[first:(first + len - (c > rem))]
end

@inline function _threaded_foreach(f, ::Union{Bool, BaseThreads}, range)
    nchunks = _nchunks(range)
    return Threads.@threads :dynamic for c in 1:nchunks
        for i in _chunk(range, c, nchunks)
            f(i, c)
        end
    end
end

function _polyester_foreach(args...)
    throw(
        ArgumentError(
            "PolyesterThreads() requires Polyester.jl to be loaded. Add `using Polyester` to your code."
        )
    )
end

@inline function _threaded_foreach(f, ::PolyesterThreads, range)
    nchunks = _nchunks(range)
    return _polyester_foreach(1:nchunks) do c
        for i in _chunk(range, c, nchunks)
            f(i, c)
        end
    end
end

macro threaded(option, ex)
    return _threaded_loop(option, gensym(:slot), ex)
end

macro threaded(option, slot, ex)
    return _threaded_loop(option, slot, ex)
end

function _threaded_loop(option, slot, ex)
    ex.head === :for || error("@threaded expects a for loop")
    loop_var = esc(ex.args[1].args[1])
    range_expr = esc(ex.args[1].args[2])
    body = esc(ex.args[2])
    slot = esc(slot)
    return quote
        option = $(esc(option))
        if isthreaded(option)
            _threaded_foreach(option, $range_expr) do $loop_var, $slot
                $body
            end
        else
            let $slot = 1
                $(esc(ex))
            end
        end
    end
end

@inline function _thread_storage_size()
    return Threads.nthreads(:default) + Threads.nthreads(:interactive)
end
