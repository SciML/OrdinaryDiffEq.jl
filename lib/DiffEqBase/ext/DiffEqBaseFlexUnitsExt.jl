module DiffEqBaseFlexUnitsExt

using DiffEqBase
import SciMLBase: unitfulvalue, value
using FlexUnits

# Support adaptive errors should be errorless for exponentiation
value(::Type{Quantity{T, U}}) where {T, U} = T
value(x::Quantity{T, U}) where {T, U} = dstrip(x)

unitfulvalue(::Type{T}) where {T <: Quantity} = T
unitfulvalue(x::Quantity) = x

DiffEqBase.stripunits(x::Quantity) = dstrip(x)

@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:Quantity,
            N,
        },
        t
    ) where {N}
    return sqrt(sum(x -> abs2(value(x)), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Quantity, N},
        t
    ) where {N}
    return sqrt(sum(x -> abs2(value(x)), u) / max(length(u), 1))
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Quantity, t) = abs(value(u))
@inline function DiffEqBase.UNITLESS_ABS2(x::Quantity)
    return real(abs2(dstrip(x)))
end

DiffEqBase._rate_prototype(u, t::Quantity, onet) = u / unit(t)
end
