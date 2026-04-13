using OrdinaryDiffEq

## https://github.com/SciML/DifferentialEquations.jl/issues/1013

mutable struct SomeObject
    position::Any
    velocity::Any
    trajectory::Any
end

object = SomeObject(0, 1, nothing)

# Current dynamics don't involve the object for the sake of MWE, but they could.
function dynamics(du, u, p, t)
    du[1] = u[2]
    return du[2] = -u[2]
end

for i in 1:2
    initial_state = [0, 0]
    tspan = (0.0, 5.0)
    prob = ODEProblem(dynamics, initial_state, tspan, object)
    sol = solve(prob, Tsit5())
    object.trajectory = sol
end

# https://github.com/SciML/DiffEqBase.jl/issues/1003

f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)

prob2 = ODEProblem((du, u, p, t) -> du[1] = 1, [0.0], (0, 10), (; x = sol))
solve(prob2, Tsit5())
