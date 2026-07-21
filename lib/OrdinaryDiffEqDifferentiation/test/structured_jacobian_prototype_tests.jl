using LinearAlgebra: Bidiagonal, Diagonal, Tridiagonal
using OrdinaryDiffEqRosenbrock: Rosenbrock23
using SciMLBase: ODEFunction, ODEProblem, init, solve
using Test

function structured_f!(du, u, p, t)
    return du .= u ./ t
end

function structured_jac!(J, u, p, t)
    J[1, 1] = 1 / t
    J[2, 2] = 1 / t
    J[1, 2] = 0
    J[2, 1] = 0
    return J
end

prototypes = (
    Diagonal(zeros(2)),
    Bidiagonal(zeros(2), zeros(1), :U),
    Bidiagonal(zeros(2), zeros(1), :L),
    Tridiagonal(zeros(1), zeros(2), zeros(1)),
)

for jac_prototype in prototypes
    f = ODEFunction(structured_f!; jac = structured_jac!, jac_prototype)
    prob = ODEProblem(f, ones(2), (1.0, 10.0))
    @test init(prob, Rosenbrock23()) !== nothing
end

f = ODEFunction(structured_f!; jac = structured_jac!, jac_prototype = first(prototypes))
prob = ODEProblem(f, ones(2), (1.0, 10.0))
sol = solve(prob, Rosenbrock23())

@test sol.u[end] ≈ [10.0, 10.0]
@test length(sol.t) < 60
