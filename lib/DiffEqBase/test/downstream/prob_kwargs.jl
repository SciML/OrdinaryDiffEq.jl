using OrdinaryDiffEq, Test
function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan, alg = Tsit5())
@test_nowarn sol = solve(prob, reltol = 1.0e-6)
sol = solve(prob, reltol = 1.0e-6)
@test sol.alg isa Tsit5

new_u0 = rand(3)
sol = solve(prob, u0 = new_u0)
@test sol.prob.u0 === new_u0
