using OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqSDIRK
using SciMLBase
using SciMLOperators
using SparseArrays
using LinearAlgebra
using Test

# ScalarOperator (λ·I) reports axes(mm) == (), unlike UniformScaling, so it fell
# through the `mass_matrix isa UniformScaling` special case and hit the
# `axes(mass_matrix) == axes(W)` boundscheck meant for full mass matrices.
J = [1.0 2.0; 3.0 4.0]
λ = 2.0
W_expected = J - λ * inv(0.5) * I

# `Matrix` specialization
W = similar(J)
OrdinaryDiffEqDifferentiation.jacobian2W!(W, ScalarOperator(λ), 0.5, J)
@test W ≈ W_expected

W_uniform = similar(J)
OrdinaryDiffEqDifferentiation.jacobian2W!(W_uniform, λ * I, 0.5, J)
@test W ≈ W_uniform

# generic `AbstractMatrix` method, reached through a sparse J
Js = sparse(J)
W_sparse = similar(Js)
W_sparse_uniform = similar(Js)
OrdinaryDiffEqDifferentiation.jacobian2W!(W_sparse, ScalarOperator(λ), 0.5, Js)
OrdinaryDiffEqDifferentiation.jacobian2W!(W_sparse_uniform, λ * I, 0.5, Js)
@test W_sparse == W_sparse_uniform
@test W_sparse ≈ W_expected

# out-of-place
@test OrdinaryDiffEqDifferentiation.jacobian2W(ScalarOperator(λ), 0.5, J) ≈ W_expected

# end-to-end: the `solve` reported in #3915
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return nothing
end
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 1.0)
sol_op = solve(
    ODEProblem(ODEFunction(lorenz!; mass_matrix = ScalarOperator(λ)), u0, tspan),
    TRBDF2(), abstol = 1.0e-10, reltol = 1.0e-10
)
sol_uniform = solve(
    ODEProblem(ODEFunction(lorenz!; mass_matrix = λ * I), u0, tspan),
    TRBDF2(), abstol = 1.0e-10, reltol = 1.0e-10
)
@test SciMLBase.successful_retcode(sol_op)
@test sol_op.u[end] ≈ sol_uniform.u[end]

# The case `UniformScaling` cannot express: a time-varying scalar mass matrix.
# M(t) u' = -u with M(t) = (1+t) I has the exact solution u(t) = u0 / (1+t), so
# a stale or ignored `update_func` shows up directly in the answer.
mm_t = ScalarOperator(1.0; update_func = (a, u, p, t) -> 1.0 + t)
decay!(du, u, p, t) = (du .= -u; nothing)
sol_t = solve(
    ODEProblem(ODEFunction(decay!; mass_matrix = mm_t), [1.0], (0.0, 1.0)),
    TRBDF2(), abstol = 1.0e-10, reltol = 1.0e-10
)
@test SciMLBase.successful_retcode(sol_t)
@test sol_t.u[end][1] ≈ 0.5 atol = 1.0e-5
