using OrdinaryDiffEq, Test
using DiffEqBase
function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
@test_nowarn sol = solve(prob, Tsit5(), reltol = 1.0e-6)
sol = solve(prob, Tsit5(), rel_tol = 1.0e-6, kwargshandle = SciMLBase.KeywordArgWarn)
@test_logs (:warn, SciMLBase.KWARGWARN_MESSAGE) sol = solve(
    prob, Tsit5(), rel_tol = 1.0e-6, kwargshandle = SciMLBase.KeywordArgWarn
)
@test_throws SciMLBase.CommonKwargError solve(prob, Tsit5(), rel_tol = 1.0e-6)

prob = ODEProblem(lorenz, u0, tspan, test = 2.0, kwargshandle = SciMLBase.KeywordArgWarn)
@test_logs (:warn, SciMLBase.KWARGWARN_MESSAGE) sol = solve(prob, Tsit5(), reltol = 1.0e-6)
