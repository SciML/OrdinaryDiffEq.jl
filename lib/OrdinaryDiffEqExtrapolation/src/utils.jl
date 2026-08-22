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

@inline function _threaded_foreach(f, ::Union{Bool, BaseThreads}, range)
    return Threads.@threads :static for i in range
        f(i)
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
    return _polyester_foreach(f, range)
end

macro threaded(option, ex)
    ex.head === :for || error("@threaded expects a for loop")
    loop_var = esc(ex.args[1].args[1])
    range_expr = esc(ex.args[1].args[2])
    body = esc(ex.args[2])
    return quote
        option = $(esc(option))
        if isthreaded(option)
            _threaded_foreach(option, $range_expr) do $loop_var
                $body
            end
        else
            $(esc(ex))
        end
    end
end

@inline function _thread_storage_size()
    return Threads.threadpoolsize(:default) + Threads.threadpoolsize(:interactive)
end
