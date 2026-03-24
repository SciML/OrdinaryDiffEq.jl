function abs2_and_sum(x, y)
    return reduce(+, x, init = zero(real(value(eltype(x))))) +
        reduce(+, y, init = zero(real(value(eltype(y)))))
end
UNITLESS_ABS2(x::Number) = abs2(x)

function UNITLESS_ABS2(x::AbstractArray{<:Number})
    return mapreduce(UNITLESS_ABS2, +, x, init = zero(real(value(eltype(x)))))
end

function UNITLESS_ABS2(x::AbstractArray)
    return mapreduce(UNITLESS_ABS2, abs2_and_sum, x, init = zero(real(value(eltype(x)))))
end
function UNITLESS_ABS2(x::RecursiveArrayTools.AbstractVectorOfArray)
    return mapreduce(UNITLESS_ABS2, abs2_and_sum, x.u, init = zero(real(value(eltype(x)))))
end
function UNITLESS_ABS2(x::RecursiveArrayTools.ArrayPartition)
    return mapreduce(UNITLESS_ABS2, abs2_and_sum, x.x, init = zero(real(value(eltype(x)))))
end

UNITLESS_ABS2(f::F, x::Number) where {F} = abs2(f(x))
function UNITLESS_ABS2(f::F, x::AbstractArray) where {F}
    return mapreduce(
        UNITLESS_ABS2 ∘ f, abs2_and_sum, x;
        init = zero(real(value(eltype(x))))
    )
end
function UNITLESS_ABS2(f::F, x::RecursiveArrayTools.ArrayPartition) where {F}
    return mapreduce(
        UNITLESS_ABS2 ∘ f, abs2_and_sum, x.x;
        init = zero(real(value(eltype(x))))
    )
end

recursive_length(u::AbstractArray{<:Number}) = length(u)
recursive_length(u::Number) = length(u)
recursive_length(u::AbstractArray{<:AbstractArray}) = sum(recursive_length, u)
recursive_length(u::RecursiveArrayTools.ArrayPartition) = sum(recursive_length, u.x)
recursive_length(u::RecursiveArrayTools.VectorOfArray) = sum(recursive_length, u.u)
function recursive_length(
        u::AbstractArray{
            <:StaticArraysCore.StaticArray{S, <:Number},
        }
    ) where {S}
    return prod(Size(eltype(u))) * length(u)
end

ODE_DEFAULT_NORM(u::Union{AbstractFloat, Complex}, t) = @fastmath abs(u)

function ODE_DEFAULT_NORM(f::F, u::Union{AbstractFloat, Complex}, t) where {F}
    return @fastmath abs(f(u))
end

function ODE_DEFAULT_NORM(u::Array{T}, t) where {T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    return Base.FastMath.sqrt_fast(real(x) / max(length(u), 1))
end

function ODE_DEFAULT_NORM(
        f::F,
        u::Union{Array{T}, Iterators.Zip{<:Tuple{Vararg{Array{T}}}}},
        t
    ) where {F, T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(f(ui))
    end
    return Base.FastMath.sqrt_fast(real(x) / max(length(u), 1))
end

function ODE_DEFAULT_NORM(
        u::StaticArraysCore.StaticArray{<:Tuple, T},
        t
    ) where {T <: Union{AbstractFloat, Complex}}
    return Base.FastMath.sqrt_fast(real(sum(abs2, u)) / max(length(u), 1))
end

function ODE_DEFAULT_NORM(
        f::F, u::StaticArraysCore.StaticArray{<:Tuple, T},
        t
    ) where {F, T <: Union{AbstractFloat, Complex}}
    return Base.FastMath.sqrt_fast(real(sum(abs2 ∘ f, u)) / max(length(u), 1))
end

function ODE_DEFAULT_NORM(
        u::Union{
            AbstractArray,
            RecursiveArrayTools.AbstractVectorOfArray,
        },
        t
    )
    return Base.FastMath.sqrt_fast(UNITLESS_ABS2(u) / max(recursive_length(u), 1))
end

function ODE_DEFAULT_NORM(f::F, u::AbstractArray, t) where {F}
    return Base.FastMath.sqrt_fast(UNITLESS_ABS2(f, u) / max(recursive_length(u), 1))
end

ODE_DEFAULT_NORM(u, t) = norm(u)
ODE_DEFAULT_NORM(f::F, u, t) where {F} = norm(f.(u))

ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false
function ODE_DEFAULT_PROG_MESSAGE(dt, u::Array, p, t)
    tmp = u[1]
    for i in eachindex(u)
        tmp = ifelse(abs(u[i]) > abs(tmp), u[i], tmp)
    end
    return "dt=" * string(dt) * "\nt=" * string(t) * "\nmax u=" * string(tmp)
end
function ODE_DEFAULT_PROG_MESSAGE(dt, u, p, t)
    return "dt=" * string(dt) * "\nt=" * string(t) * "\nmax u=" * string(maximum(abs.(u)))
end

NAN_CHECK(x::Number) = isnan(x)
NAN_CHECK(x::Enum) = false
function NAN_CHECK(x::Union{AbstractArray, RecursiveArrayTools.AbstractVectorOfArray})
    return any(
        NAN_CHECK, x
    )
end
NAN_CHECK(x::RecursiveArrayTools.ArrayPartition) = any(NAN_CHECK, x.x)

INFINITE_OR_GIANT(x::Number) = !isfinite(x)
function INFINITE_OR_GIANT(
        x::Union{
            AbstractArray, RecursiveArrayTools.AbstractVectorOfArray,
        }
    )
    return any(
        INFINITE_OR_GIANT, x
    )
end
INFINITE_OR_GIANT(x::RecursiveArrayTools.ArrayPartition) = any(INFINITE_OR_GIANT, x.x)
ODE_DEFAULT_UNSTABLE_CHECK(dt, u, p, t) = false
function ODE_DEFAULT_UNSTABLE_CHECK(dt, u::Union{Number, AbstractArray{<:Number}}, p, t)
    return INFINITE_OR_GIANT(u)
end
