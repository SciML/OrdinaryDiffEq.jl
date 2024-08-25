using OrdinaryDiffEq, Test
import SciMLBase

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 0.5)
prob = ODEProblem(lorenz!, u0, tspan)

sol = solve(prob, Rosenbrock23())

@test isnothing(SciMLBase.strip_interpolation(sol.interp).f)
@test isnothing(SciMLBase.strip_interpolation(sol.interp).cache.jac_config)
