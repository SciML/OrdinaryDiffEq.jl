using OrdinaryDiffEq, Test, ADTypes

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
