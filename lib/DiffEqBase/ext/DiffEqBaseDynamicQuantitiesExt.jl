module DiffEqBaseDynamicQuantitiesExt

using DiffEqBase
using DynamicQuantities
using LinearAlgebra
import DiffEqBase: default_factorize
import SciMLBase: value, unitfulvalue
import RecursiveArrayTools: recursive_unitless_bottom_eltype, recursive_unitless_eltype

@inline DiffEqBase.ODE_DEFAULT_NORM(u::UnionAbstractQuantity, t) = abs(ustrip(u))
@inline function DiffEqBase.UNITLESS_ABS2(x::UnionAbstractQuantity)
    return real(abs2(ustrip(x)))
end

DiffEqBase.stripunits(x::UnionAbstractQuantity) = ustrip(x)

# `value` strips everything (units + AD/uncertainty) by recursing on the stripped
# numeric value; `unitfulvalue` strips AD but keeps units. Parallels the Unitful
# and FlexUnits extensions. Without these, `DiffEqBase.value(t::Quantity)` was the
# identity default from SciMLBase and `eps(value(t))` in initdt.jl hit
# `MethodError: no method matching eps(::Quantity)`.
value(x::Type{<:UnionAbstractQuantity}) = value(DynamicQuantities.value_type(x))
value(x::UnionAbstractQuantity) = value(ustrip(x))
unitfulvalue(x::Type{T}) where {T <: UnionAbstractQuantity} = T
unitfulvalue(x::UnionAbstractQuantity) = x

# RecursiveArrayTools' generic `Number` dispatch uses `typeof(one(T))`, but DynamicQuantities'
# `one(::Type{<:Quantity})` returns a dimensionless Quantity, so units survive. Recurse into
# the value type instead so downstream `zero(uBottomEltypeNoUnits)` lands on a bare scalar.
recursive_unitless_bottom_eltype(::Type{T}) where {T <: UnionAbstractQuantity} =
    recursive_unitless_bottom_eltype(DynamicQuantities.value_type(T))
recursive_unitless_eltype(::Type{T}) where {T <: UnionAbstractQuantity} =
    recursive_unitless_eltype(DynamicQuantities.value_type(T))

DiffEqBase._rate_prototype(u, t::UnionAbstractQuantity, onet) = u / oneunit(t)
DiffEqBase.timedepentdtmin(t::UnionAbstractQuantity, dtmin) =
    abs(ustrip(dtmin / oneunit(t)) * oneunit(t))

# Rosenbrock/SDIRK solvers form W/J matrices with Quantity eltype. Factorize/solve in
# value-space (Float64), but return solutions with the RHS units.
struct DQUnitlessLU{F, UT}
    F::F
    ut::UT
end

@inline function _infer_ut(A::AbstractMatrix{<:UnionAbstractQuantity})
    @inbounds for a in A
        va = ustrip(a)
        if !iszero(va)
            return oneunit(inv(a))
        end
    end
    return oneunit(1.0)
end

function default_factorize(A::AbstractMatrix{<:UnionAbstractQuantity})
    isempty(A) && return DQUnitlessLU(
        lu(Matrix{Float64}(undef, 0, 0); check = false),
        oneunit(1.0),
    )
    ut = _infer_ut(A)
    return DQUnitlessLU(lu(ustrip.(A); check = false), ut)
end

function LinearAlgebra.ldiv!(
        x::AbstractVector{<:UnionAbstractQuantity},
        W::DQUnitlessLU,
        b::AbstractVector{<:UnionAbstractQuantity},
    )
    vb = ustrip.(b)
    vx = similar(vb)
    LinearAlgebra.ldiv!(vx, W.F, vb)
    @inbounds for i in eachindex(x)
        x[i] = vx[i] * (oneunit(b[i]) * W.ut)
    end
    return x
end

function Base.:(\)(W::DQUnitlessLU, b::AbstractVector{<:UnionAbstractQuantity})
    vb = ustrip.(b)
    vx = W.F \ vb
    out = similar(b)
    @inbounds for i in eachindex(out)
        out[i] = vx[i] * (oneunit(b[i]) * W.ut)
    end
    return out
end

end
