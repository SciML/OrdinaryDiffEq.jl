using OrdinaryDiffEqDifferentiation
using SciMLOperators
using LinearAlgebra
using Test

# ScalarOperator (λ·I) reports axes(mm) == (), unlike UniformScaling, so it fell
# through the `mass_matrix isa UniformScaling` special case and hit the
# `axes(mass_matrix) == axes(W)` boundscheck meant for full mass matrices.
J = [1.0 2.0; 3.0 4.0]
λ = 2.0
W_expected = J - λ * inv(0.5) * I

W = similar(J)
OrdinaryDiffEqDifferentiation.jacobian2W!(W, ScalarOperator(λ), 0.5, J)
@test W ≈ W_expected

W_uniform = similar(J)
OrdinaryDiffEqDifferentiation.jacobian2W!(W_uniform, λ * I, 0.5, J)
@test W ≈ W_uniform

@test OrdinaryDiffEqDifferentiation.jacobian2W(ScalarOperator(λ), 0.5, J) ≈ W_expected

W_dense = Matrix(J)
OrdinaryDiffEqDifferentiation.jacobian2W!(W_dense, ScalarOperator(λ), 0.5, Matrix(J))
@test W_dense ≈ W_expected
