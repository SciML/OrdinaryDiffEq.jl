"""
$(TYPEDEF)

Holds a tableau which defines an explicit Runge-Kutta method.
"""
mutable struct ExplicitRKTableau{MType <: AbstractMatrix, VType <: AbstractVector,
        CType <: AbstractVector, S, IType} <:
    ODERKTableau
    A::MType
    c::CType
    α::VType
    αEEst::VType
    d::VType # dense output coefficients
    stages::Int
    order::Int
    adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
    fsal::Bool
    stability_size::S
    B_interp::IType
end
function ExplicitRKTableau(
        A::MType, c::CType, α::VType, order;
        adaptiveorder = 0, αEEst = similar(α, 0),
        fsal = false, stability_size = 0.0,
        d = similar(α, 0), B_interp::IType = nothing
    ) where {MType, VType, CType, IType}
    S = typeof(stability_size)
    return ExplicitRKTableau{MType, VType, CType, S, IType}(
        A, c, α, αEEst, d, length(α), order, adaptiveorder,
        fsal, stability_size, B_interp
    )
end

"""
$(TYPEDEF)

Holds a tableau which defines an implicit Runge-Kutta method.
"""
mutable struct ImplicitRKTableau{MType <: AbstractMatrix, VType <: AbstractVector,
        CType <: AbstractVector} <:
    ODERKTableau
    A::MType
    c::CType
    α::VType
    αEEst::VType
    stages::Int
    order::Int
    adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
end
function ImplicitRKTableau(
        A::MType, c::CType, α::VType, order;
        adaptiveorder = 0, αEEst = similar(α, 0)
    ) where {MType, VType, CType}
    return ImplicitRKTableau{MType, VType, CType}(
        A, c, α, αEEst, length(α), order, adaptiveorder)
end
