using LinearAlgebra
using SparseArrays
using CUDA
using CUDA.CUSPARSE
using OrdinaryDiffEq

# Linear exponential solvers
A = MatrixOperator([2.0 -1.0; -1.0 2.0])
u0 = ones(2)

A_gpu = MatrixOperator(cu([2.0 -1.0; -1.0 2.0]))
u0_gpu = cu(ones(2))
prob_gpu = ODEProblem(A_gpu, u0_gpu, (0.0, 1.0))

sol_analytic = exp(1.0 * Matrix(A)) * u0

@test_broken sol1_gpu = solve(prob_gpu, LinearExponential(krylov = :off))(1.0) |> Vector
sol2_gpu = solve(prob_gpu, LinearExponential(krylov = :simple))(1.0) |> Vector
sol3_gpu = solve(prob_gpu, LinearExponential(krylov = :adaptive))(1.0) |> Vector

@test_broken isapprox(sol1_gpu, sol_analytic, rtol = 1e-6)
@test isapprox(sol2_gpu, sol_analytic, rtol = 1e-6)
@test isapprox(sol3_gpu, sol_analytic, rtol = 1e-6)

A2_gpu = MatrixOperator(cu(sparse([2.0 -1.0; -1.0 2.0])))
prob2_gpu = ODEProblem(A2_gpu, u0_gpu, (0.0, 1.0))

@test_broken sol2_1_gpu = solve(prob2_gpu, LinearExponential(krylov = :off))(1.0) |> Vector
sol2_2_gpu = solve(prob2_gpu, LinearExponential(krylov = :simple))(1.0) |> Vector
sol2_3_gpu = solve(prob2_gpu, LinearExponential(krylov = :adaptive))(1.0) |> Vector

@test_broken isapprox(sol2_1_gpu, sol_analytic, rtol = 1e-6)
@test isapprox(sol2_2_gpu, sol_analytic, rtol = 1e-6)
@test isapprox(sol2_3_gpu, sol_analytic, rtol = 1e-6)
