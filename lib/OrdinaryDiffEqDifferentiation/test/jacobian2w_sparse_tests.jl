using OrdinaryDiffEqDifferentiation
using LinearAlgebra
using SparseArrays
using Test

# `jacobian2W!` forms W = J - M/dtgamma. With a scalar mass matrix it writes the
# diagonal, which a sparse GPU `W` supports neither by scalar indexing nor by
# broadcasting into a diagonal view, so that case takes the allocating path instead.
# The allocating path is only selected for GPU storage, so the equivalence of the two
# formulas is what can be checked on CPU; the GPU dispatch itself is covered by
# test/gpu/sparse_jac_prototype.jl.

J = sparse([3.0 1.0; 0.5 2.0])
dtgamma = 0.25
invdtgamma = inv(dtgamma)

for mm in (I, 2.0 * I)
    λ = -OrdinaryDiffEqDifferentiation._scalar_massmatrix_λ(mm)

    W = similar(J)
    fill!(nonzeros(W), 0)
    OrdinaryDiffEqDifferentiation.jacobian2W!(W, mm, dtgamma, J)

    # The formula the sparse-GPU branch uses must agree with the diagonal write.
    allocating = J + (λ * invdtgamma) * I
    @test Matrix(W) ≈ Matrix(allocating)
    @test Matrix(W) ≈ Matrix(J) - OrdinaryDiffEqDifferentiation._scalar_massmatrix_λ(mm) *
        invdtgamma * Matrix(I, 2, 2)
end

W = copy(J)
@test OrdinaryDiffEqDifferentiation._update_sparse_diagonal!(
    W, -1.0, invdtgamma, J
)
@test Matrix(W) ≈ Matrix(J) - invdtgamma * Matrix(I, 2, 2)

Jmissing = sparse([1, 2], [2, 1], [2.0, 3.0], 2, 2)
Wmissing = copy(Jmissing)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wmissing, I, dtgamma, Jmissing)
@test Matrix(Wmissing) ≈ Matrix(Jmissing) - invdtgamma * Matrix(I, 2, 2)
@test nnz(Wmissing) == nnz(Jmissing) + 2

# Partial diagonal: some columns take the fast write, others the insertion fallback.
Jpartial = sparse([1, 1, 2], [1, 2, 1], [3.0, 1.0, 0.5], 2, 2)
Wpartial = copy(Jpartial)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wpartial, I, dtgamma, Jpartial)
@test Matrix(Wpartial) ≈ Matrix(Jpartial) - invdtgamma * Matrix(I, 2, 2)
@test nnz(Wpartial) == nnz(Jpartial) + 1

# CPU sparse storage is scalar-indexable, so it keeps the in-place diagonal write; a
# dense matrix is not sparse at all. Only GPU storage takes the allocating branch.
W = similar(J)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(W)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(ones(2, 2))

# The dense `Matrix` method must be unaffected by the sparse guard.
Jd = Matrix(J)
Wd = similar(Jd)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wd, I, dtgamma, Jd)
@test Wd ≈ Jd - invdtgamma * Matrix(I, 2, 2)
