using Test
using OrdinaryDiffEq
using LinearAlgebra, LinearSolve

import OrdinaryDiffEqCore.dolinsolve

n = 8
dt = 1 / 1000
u0 = ones(n)
tspan = (0.0, 1.0)

M1 = 2ones(n) |> Diagonal #|> Array
M2 = 2ones(n) |> Diagonal #|> Array

f1 = (du, u, p, t) -> du .= M1 * u
f2 = (du, u, p, t) -> du .= M2 * u
prob = SplitODEProblem(f1, f2, u0, tspan)

for algname in (
        :SBDF2,
        :SBDF3,
        :KenCarp47,
    )
    @testset "$algname" begin
        alg0 = @eval $algname()
        alg1 = @eval $algname(linsolve = LUFactorization())

        kwargs = (dt = dt,)

        solve(prob, alg0; kwargs...)
        @test SciMLBase.__solve(prob, alg0; kwargs...).retcode == ReturnCode.Success
        @test SciMLBase.__solve(prob, alg1; kwargs...).retcode == ReturnCode.Success
    end
end

f1 = M1 |> MatrixOperator
f2 = M2 |> MatrixOperator
prob = SplitODEProblem(f1, f2, u0, tspan)

for algname in (
        :SBDF2,
        :SBDF3,
        :KenCarp47,
    )
    @testset "$algname" begin
        alg0 = @eval $algname()

        kwargs = (dt = dt,)

        solve(prob, alg0; kwargs...)
        @test SciMLBase.__solve(prob, alg0; kwargs...).retcode == ReturnCode.Success
    end
end

###
# custom linsolve function
###

function linsolve(A, b, u, p, newA, Pl, Pr, solverdata; kwargs...)
    # avoid preconditioner monkeybusiness
    prob = LinearProblem(A, b; u0 = u)
    solve(prob, nothing)
    return u
end

alg = KenCarp47(linsolve = LinearSolveFunction(linsolve))

@test solve(prob, alg).retcode == ReturnCode.Success
