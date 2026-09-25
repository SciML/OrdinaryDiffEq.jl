using OrdinaryDiffEqDifferentiation
using SparseArrays
using Test

# CPU sparse path uses broadcast; verify W = J_u + cj * J_du numerically.
# GPU CSR (allocating path) is covered by monorepo test/gpu/dae_tests.jl.
J_u = sparse([1.0 0.0; 0.5 2.0])
J_du = sparse([0.0 1.0; 0.0 0.5])
cj = 3.0
W = similar(J_u)
fill!(nonzeros(W), 0)
OrdinaryDiffEqDifferentiation.dae_jacobian2W!(W, J_u, J_du, cj)
@test Matrix(W) ≈ Matrix(J_u) + cj * Matrix(J_du)

W_oop = OrdinaryDiffEqDifferentiation.dae_jacobian2W(J_u, J_du, cj)
@test Matrix(W_oop) ≈ Matrix(J_u) + cj * Matrix(J_du)

# CPU sparse takes the broadcast path (fast_scalar_indexing storage); a dense
# matrix is not sparse at all. Real cuSPARSE types subtype AbstractSparseMatrix
# with GPU nonzeros storage and take the allocating path (covered on GPU CI).
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(W)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(ones(2, 2))
