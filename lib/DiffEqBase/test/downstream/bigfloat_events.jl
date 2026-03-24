using OrdinaryDiffEq

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = BigFloat[1.0; 0.0; 0.0]
tspan = (big(0.0), big(100.0))
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5(), save_everystep = false)

x = sol.u[end]

import LinearAlgebra.norm
function condition(u, t, integrator)
    return norm(u - x) - 0.1
end

affect!(integrator) = terminate!(integrator)

cb = ContinuousCallback(condition, affect!)

sol2 = solve(prob, Tsit5(), save_everystep = false, callback = cb)
