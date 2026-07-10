using OrdinaryDiffEqDifferentiation
using SparseArrays
using Test

# CPU sparse path uses broadcast; verify W = J_u + cj * J_du numerically.
# GPU CSR is covered by monorepo test/gpu/dae_tests.jl (uses nzVal detection).
J_u = sparse([1.0 0.0; 0.5 2.0])
J_du = sparse([0.0 1.0; 0.0 0.5])
cj = 3.0
W = similar(J_u)
fill!(nonzeros(W), 0)
OrdinaryDiffEqDifferentiation.dae_jacobian2W!(W, J_u, J_du, cj)
@test Matrix(W) ≈ Matrix(J_u) + cj * Matrix(J_du)

W_oop = OrdinaryDiffEqDifferentiation.dae_jacobian2W(J_u, J_du, cj)
@test Matrix(W_oop) ≈ Matrix(J_u) + cj * Matrix(J_du)

# Duck-type path for types that look like cuSPARSE (have .nzVal) but are not
# AbstractSparseMatrix — mirrors CuSparseMatrixCSR detection.
struct FakeCuSparse{T} <: AbstractMatrix{T}
    data::Matrix{T}
    nzVal::Vector{T}
end
Base.size(A::FakeCuSparse) = size(A.data)
Base.getindex(A::FakeCuSparse, i::Int, j::Int) = A.data[i, j]
Base.setindex!(A::FakeCuSparse, v, i::Int, j::Int) = (A.data[i, j] = v)
function Base.copyto!(dest::FakeCuSparse, src::AbstractMatrix)
    dest.data .= src
    return dest
end
Base.:+(A::FakeCuSparse, B::AbstractMatrix) = FakeCuSparse(A.data + Matrix(B), vec(A.data + Matrix(B)))
Base.:+(A::AbstractMatrix, B::FakeCuSparse) = FakeCuSparse(Matrix(A) + B.data, vec(Matrix(A) + B.data))
Base.:+(A::FakeCuSparse, B::FakeCuSparse) = FakeCuSparse(A.data + B.data, vec(A.data + B.data))
Base.:*(a::Number, A::FakeCuSparse) = FakeCuSparse(a * A.data, a * A.nzVal)

@test OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(
    FakeCuSparse(ones(2, 2), ones(4))
)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(ones(2, 2))

J_u_f = FakeCuSparse([1.0 0.0; 0.5 2.0], [1.0, 0.5, 2.0])
J_du_f = FakeCuSparse([0.0 1.0; 0.0 0.5], [1.0, 0.5])
W_f = FakeCuSparse(zeros(2, 2), zeros(4))
OrdinaryDiffEqDifferentiation.dae_jacobian2W!(W_f, J_u_f, J_du_f, cj)
@test W_f.data ≈ [1.0 3.0; 0.5 3.5]
