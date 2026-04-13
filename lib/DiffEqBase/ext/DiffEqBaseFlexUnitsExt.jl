module DiffEqBaseFlexUnitsExt

using DiffEqBase
import SciMLBase: unitfulvalue, value
using FlexUnits

# Support adaptive errors should be errorless for exponentiation
value(::Type{Quantity{T, U}}) where {T, U} = T
value(x::Quantity{T, U}) where {T, U} = dstrip(x)

unitfulvalue(::Type{T}) where {T <: Quantity} = T
unitfulvalue(x::Quantity) = x

@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:Quantity,
            N,
        },
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Quantity, N},
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Quantity, t) = abs(value(u))
@inline function DiffEqBase.UNITLESS_ABS2(x::Quantity)
    return real(abs2(dstrip(x)))
end

DiffEqBase._rate_prototype(u, t::Quantity, onet) = u / unit(t)
end
