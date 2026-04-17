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
# OrdinaryDiffEqLinearExponentialAlgorithm subtypes (MagnusGL6, etc.)
# have no `autodiff` field; _alg_autodiff and prepare_alg must not crash.
using OrdinaryDiffEqDifferentiation: _alg_autodiff
using OrdinaryDiffEqCore: OrdinaryDiffEqLinearExponentialAlgorithm
using DiffEqBase: prepare_alg

struct MockMagnusAlg <: OrdinaryDiffEqLinearExponentialAlgorithm
    krylov::Bool
    m::Int
    iop::Int
end

@testset "LinearExponentialAlgorithm autodiff traits (issue #3232)" begin
    mock = MockMagnusAlg(false, 30, 0)

    # _alg_autodiff must return Val{false}() instead of accessing alg.autodiff
    @test _alg_autodiff(mock) == Val{false}()

    # prepare_alg must return the algorithm unchanged (no AD preparation needed)
    u0 = ones(2)
    mock_prob = ODEProblem((du, u, p, t) -> du .= 0, u0, (0.0, 1.0))
    @test prepare_alg(mock, u0, nothing, mock_prob) === mock

    # forwarddiffs_model must return false
    @test SciMLBase.forwarddiffs_model(mock) == false
end
