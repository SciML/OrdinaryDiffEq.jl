"""
    ComplexShiftOperator{T} <: AbstractSciMLOperator{T}

Lazy operator representing `J - c * M`, where `J` is a Jacobian (or a matrix-free
Jacobian-vector-product operator such as `JVPCache`), `M` is the mass matrix and `c` is a
complex shift.

FIRK methods decouple the stage system through the W-transform into one real-shifted and
`(s - 1) ÷ 2` complex-shifted linear systems. The real-shifted system is handled by
`SciMLOperators.WOperator`; this operator covers the complex-shifted ones. Since `J` is
real, it is applied separately to the real and imaginary parts of the argument, so a real
matrix-free Jacobian can back a complex linear solve.
"""
mutable struct ComplexShiftOperator{T, MMType, CType, JType, VType} <:
    AbstractSciMLOperator{T}
    mass_matrix::MMType
    shift::CType
    J::JType
    vr::VType
    vi::VType
    Jvr::VType
    Jvi::VType
end

function ComplexShiftOperator(mass_matrix, shift, J, u)
    v = similar(_vec(u), eltype(J))
    fill!(v, zero(eltype(v)))
    T = Complex{real(eltype(J))}
    return ComplexShiftOperator{
        T, typeof(mass_matrix), typeof(shift), typeof(J), typeof(v),
    }(
        mass_matrix, shift, J, v, similar(v), similar(v), similar(v)
    )
end

Base.size(W::ComplexShiftOperator) = size(W.J)
Base.size(W::ComplexShiftOperator, d::Integer) = d <= 2 ? size(W)[d] : 1
SciMLBase.isinplace(::ComplexShiftOperator) = true
SciMLOperators.isconvertible(::ComplexShiftOperator) = false
SciMLOperators.has_ldiv(::ComplexShiftOperator) = false
SciMLOperators.has_ldiv!(::ComplexShiftOperator) = false

function SciMLOperators.update_coefficients!(
        W::ComplexShiftOperator, u = nothing, p = nothing, t = nothing; shift = nothing
    )
    if (u !== nothing) && (p !== nothing) && (t !== nothing)
        SciMLOperators.update_coefficients!(W.J, u, p, t)
        SciMLOperators.update_coefficients!(W.mass_matrix, u, p, t)
    end
    shift !== nothing && (W.shift = shift)
    return W
end

function LinearAlgebra.mul!(
        Y::AbstractVector, W::ComplexShiftOperator, B::AbstractVector
    )
    (; vr, vi, Jvr, Jvi, mass_matrix, shift) = W
    @.. vr = real(B)
    @.. vi = imag(B)
    mul!(Jvr, W.J, vr)
    mul!(Jvi, W.J, vi)
    if mass_matrix isa UniformScaling
        c = shift * mass_matrix.λ
        @.. Y = Jvr + im * Jvi - c * B
    else
        mul!(Y, mass_matrix, B)
        @.. Y = Jvr + im * Jvi - shift * Y
    end
    return Y
end

matrixfree_jacobian(W::WOperator) = W.jacvec === nothing ? W.J : W.jacvec

function firk_update_W!(
        W1::WOperator, W2::ComplexShiftOperator, J, mass_matrix, γdt, αdt, βdt, uprev, p, t
    )
    update_coefficients!(W1, uprev, p, t; gamma = inv(γdt))
    update_coefficients!(W2, uprev, p, t; shift = αdt + βdt * im)
    return nothing
end

function firk_update_W!(W1, W2, J, mass_matrix, γdt, αdt, βdt, uprev, p, t)
    @inbounds for II in CartesianIndices(J)
        W1[II] = -γdt * mass_matrix[Tuple(II)...] + J[II]
        W2[II] = -(αdt + βdt * im) * mass_matrix[Tuple(II)...] + J[II]
    end
    return nothing
end

function firk_matrixfree_unsupported(algname)
    return error(
        "Non-concrete (matrix-free) Jacobians are not yet supported by $algname. " *
            "Matrix-free `linsolve` choices such as `KrylovJL_GMRES()` are currently " *
            "supported by `RadauIIA5` only; use `RadauIIA5(linsolve = KrylovJL_GMRES())` " *
            "for a matrix-free FIRK method, or pass a factorization-based `linsolve` " *
            "(the default, or e.g. `LUFactorization()`/`KLUFactorization()` with a " *
            "`jac_prototype`) to $algname."
    )
end
