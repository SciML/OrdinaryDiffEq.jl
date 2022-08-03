using Test
using OrdinaryDiffEq
using LinearAlgebra, LinearSolve
using SciMLOperators

import OrdinaryDiffEq.dolinsolve

n = 8
dt = 1 / 16
u0 = ones(n)
tspan = (0.0, 1.0)

M1 = 2ones(n) |> Diagonal #|> Array
M2 = 2ones(n) |> Diagonal #|> Array

f1 = M1 |> MatrixOperator
f2 = M2 |> MatrixOperator
prob = SplitODEProblem(f1, f2, u0, tspan)

for algname in (:SBDF2,
                :SBDF3,
                :KenCarp47)
    @testset "$algname" begin
        alg0 = @eval $algname()
        alg1 = @eval $algname(linsolve = GenericFactorization())

        kwargs = (dt = dt,)

        # expected error message
        msg = "Split ODE problem do not work with factorization linear solvers. Bug detailed in https://github.com/SciML/OrdinaryDiffEq.jl/pull/1643. Defaulting to linsolve=KrylovJL()"
        @test_logs (:warn, msg) solve(prob, alg0; kwargs...)
        @test DiffEqBase.__solve(prob, alg0; kwargs...).retcode == :Success
        @test_broken DiffEqBase.__solve(prob, alg1; kwargs...).retcode == :Success
    end
end

#####
# deep dive
#####

alg0 = KenCarp47()                                # passing case
alg1 = KenCarp47(linsolve = GenericFactorization()) # failing case

## objects
ig0 = SciMLBase.init(prob, alg0; dt = dt)
ig1 = SciMLBase.init(prob, alg1; dt = dt)

nl0 = ig0.cache.nlsolver
nl1 = ig1.cache.nlsolver

lc0 = nl0.cache.linsolve
lc1 = nl1.cache.linsolve

W0 = lc0.A
W1 = lc1.A

# perform first step
OrdinaryDiffEq.loopheader!(ig0)
OrdinaryDiffEq.loopheader!(ig1)

OrdinaryDiffEq.perform_step!(ig0, ig0.cache)
OrdinaryDiffEq.perform_step!(ig1, ig1.cache)

@test !OrdinaryDiffEq.nlsolvefail(nl0)
@test OrdinaryDiffEq.nlsolvefail(nl1)

# check operators
@test W0._concrete_form != W1._concrete_form
@test_broken W0._func_cache == W1._func_cache

# check operator application
b = ones(n)
@test W0 * b == W1 * b
@test mul!(rand(n), W0, b) == mul!(rand(n), W1, b)
#@test W0 \ b == W1 \ b

# check linear solve
lc0.b .= 1.0
lc1.b .= 1.0

solve(lc0)
solve(lc1)

@test_broken lc0.u == lc1.u

# solve contried problem using OrdinaryDiffEq machinery
linres0 = dolinsolve(ig0, lc0; A = W0, b = b, linu = ones(n), reltol = 1e-8)
linres1 = dolinsolve(ig1, lc1; A = W1, b = b, linu = ones(n), reltol = 1e-8)

@test_broken linres0 == linres1

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

@test solve(prob, alg).retcode == :Success

nothing
