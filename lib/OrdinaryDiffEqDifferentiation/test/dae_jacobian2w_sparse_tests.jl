using OrdinaryDiffEqDifferentiation
using SparseArrays
using Test

# CPU sparse path fills `nonzeros(W)` over a fixed pattern; verify
# W = J_u + cj * J_du numerically and bitwise against `muladd` broadcast.
# GPU CSR (allocating path) is covered by monorepo test/gpu/dae_tests.jl.
J_u = sparse([1.0 0.0; 0.5 2.0])
J_du = sparse([0.0 1.0; 0.0 0.5])
cj = 3.0
W = similar(J_u)
fill!(nonzeros(W), 0)
OrdinaryDiffEqDifferentiation.dae_jacobian2W!(W, J_u, J_du, cj)
@test Matrix(W) ≈ Matrix(J_u) + cj * Matrix(J_du)
@test map(bitstring, Matrix(W)) == map(bitstring, muladd.(cj, Matrix(J_du), Matrix(J_u)))

W_oop = OrdinaryDiffEqDifferentiation.dae_jacobian2W(J_u, J_du, cj)
@test Matrix(W_oop) ≈ Matrix(J_u) + cj * Matrix(J_du)

# W's pattern is widened once to `W ∪ J_u ∪ J_du` and then stays fixed, even
# when stored values become exactly 0.0 and broadcast would have dropped them.
Wstable = similar(J_u)
fill!(nonzeros(Wstable), 0)
nnz_hist = Int[]
for val in (1.0, 0.0, -2.0)
    J_du[1, 2] = val
    OrdinaryDiffEqDifferentiation.dae_jacobian2W!(Wstable, J_u, J_du, cj)
    push!(nnz_hist, nnz(Wstable))
end
@test allequal(nnz_hist)
@test nnz(Wstable) >= max(nnz(J_u), nnz(J_du))
@test Matrix(Wstable) ≈ Matrix(J_u) + cj * Matrix(J_du)

# `_sparse_pattern_stable_sub!` (the `J_du = J_du - J_u` in `calc_J_dae!`) keeps
# `J_du`'s pattern while subtracting in place.
A = sparse([1, 1, 2], [1, 2, 2], [4.0, 1.0, 0.5], 2, 2)
B = sparse([1, 2, 2], [2, 1, 2], [0.25, 1.0, 0.25], 2, 2)
Adense = Matrix(A) - Matrix(B)
@test OrdinaryDiffEqDifferentiation._sparse_pattern_stable_sub!(A, B)
@test Matrix(A) ≈ Adense
@test map(bitstring, Matrix(A)) == map(bitstring, Adense)
@test nnz(A) == 4  # A ∪ B
A2 = sparse([1, 1, 2], [1, 2, 2], [4.0, 0.0, 0.5], 2, 2)
@test OrdinaryDiffEqDifferentiation._sparse_pattern_stable_sub!(A2, B)
@test nnz(A2) == nnz(A)

# CPU sparse takes the in-place path (fast_scalar_indexing storage); a dense
# matrix is not sparse at all. Real cuSPARSE types subtype AbstractSparseMatrix
# with GPU nonzeros storage and take the allocating path (covered on GPU CI).
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(W)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(ones(2, 2))
