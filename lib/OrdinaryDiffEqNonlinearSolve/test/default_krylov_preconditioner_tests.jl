using LinearAlgebra
using OrdinaryDiffEqSDIRK
using SciMLBase
using SciMLOperators: MatrixOperator
using Test

n = 8
dx = 32.0 / n
A = zeros(n, n)
for i in 1:n
    A[i, mod1(i - 1, n)] = 10 / dx^2
    A[i, i] = -20 / dx^2
    A[i, mod1(i + 1, n)] = 10 / dx^2
end

nonlinear!(du, u, p, t) = @. du = -10u^3
u0 = cos.(2pi .* (0:(n - 1)) ./ n)
prob = SplitODEProblem(MatrixOperator(A), nonlinear!, u0, (0.0, 1.0))
integrator = init(prob, KenCarp4(); abstol = 1.0e-10, reltol = 1.0e-7)
linsolve = integrator.cache.nlsolver.cache.linsolve
weight = vec(integrator.cache.nlsolver.cache.weight)
weight .= 2:(n + 1)
left = similar(weight)
right = similar(weight)

ldiv!(left, linsolve.Pl, ones(n))
ldiv!(right, linsolve.Pr, ones(n))

@test left ≈ weight
@test right ≈ inv.(weight)

weight[1] = 0
ldiv!(left, linsolve.Pl, ones(n))
ldiv!(right, linsolve.Pr, ones(n))
expected_weight = copy(weight)
expected_weight[1] = 1
@test left ≈ expected_weight
@test right ≈ inv.(expected_weight)
