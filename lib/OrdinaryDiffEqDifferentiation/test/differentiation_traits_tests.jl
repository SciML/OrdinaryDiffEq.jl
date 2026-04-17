using OrdinaryDiffEqRosenbrock, Test, ADTypes

jac_called = Ref(false)
tgrad_called = Ref(false)

function Lotka(du, u, p, t)
    du[1] = u[1] - u[1] * u[2] # REPL[7], line 3:
    du[2] = -3 * u[2] + 1 * u[1] * u[2]
    return nothing
end

function Lotka_jac(J, u, p, t)
    jac_called.x = true
    J[1, 1] = 1.0 - u[2]
    J[1, 2] = -u[1]
    J[2, 1] = 1 * u[2]
    J[2, 2] = -3 + u[1]
    return nothing
end

function Lotka_tgrad(grad, u, p, t)
    tgrad_called.x = true
    grad[1] = 1 * 0
    grad[2] = 1 * 0
    return nothing
end

Lotka_f = ODEFunction(Lotka, jac = Lotka_jac, tgrad = Lotka_tgrad)
prob = ODEProblem(Lotka_f, ones(2), (0.0, 10.0))

good_sol = solve(prob, Rosenbrock23())

@test jac_called[]
@test tgrad_called[]

prob2 = ODEProblem(Lotka, ones(2), (0.0, 10.0))

sol = solve(prob2, Rosenbrock23(autodiff = AutoForwardDiff()))
@test ≈(good_sol[:, end], sol[:, end], rtol = 1.0e-2)

sol = solve(prob2, Rosenbrock23(autodiff = AutoForwardDiff(chunksize = 1)))
@test ≈(good_sol[:, end], sol[:, end], rtol = 1.0e-2)

sol = solve(prob2, Rosenbrock23(autodiff = AutoFiniteDiff()))
@test ≈(good_sol[:, end], sol[:, end], rtol = 1.0e-2)

# Regression test for issue #3232:
# MagnusGL6 (and all OrdinaryDiffEqLinearExponentialAlgorithm subtypes)
# have no `autodiff` field. When OrdinaryDiffEqDifferentiation is loaded,
# _alg_autodiff must not crash by trying to access alg.autodiff.
using OrdinaryDiffEqLinear
using SciMLOperators: MatrixOperator

@testset "MagnusGL6 solve with Differentiation loaded (issue #3232)" begin
    function update_func!(A, u, p, t)
        A[1, 1] = cos(t)
        A[2, 1] = sin(t)
        A[1, 2] = -sin(t)
        A[2, 2] = cos(t)
    end
    A = MatrixOperator(ones(2, 2), update_func! = update_func!)
    prob = ODEProblem(A, ones(2), (1.0, 6.0))

    # This would crash with FieldError before the fix
    sol = solve(prob, MagnusGL6(), dt = 1 / 10)
    @test sol.retcode == ReturnCode.Success
    @test length(sol.t) > 1
end
