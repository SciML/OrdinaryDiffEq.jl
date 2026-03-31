using OrdinaryDiffEq
using LabelledArrays
using ADTypes

function f(out, du, u, p, t)
    out.x = -0.04u.x + 1.0e4 * u.y * u.z - du.x
    out.y = +0.04u.x - 3.0e7 * u.y^2 - 1.0e4 * u.y * u.z - du.y
    return out.z = u.x + u.y + u.z - 1.0
end

u₀ = LVector(x = 1.0, y = 0.0, z = 0.0)
du₀ = LVector(x = -0.04, y = 0.04, z = 0.0)
tspan = (0.0, 100000.0)

differential_vars = LVector(x = true, y = true, z = false)
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars)

sol = solve(prob, DImplicitEuler())

function f1(du, u, p, t)
    du.x .= -1 .* u.x .* u.y .* p[1]
    return du.y .= -1 .* u.y .* p[2]
end
const n = 4
u_0 = @LArray fill(1000.0, 2 * n) (x = (1:n), y = ((n + 1):(2 * n)))
p = [0.1, 0.1]
prob1 = ODEProblem(f1, u_0, (0, 100.0), p)
sol = solve(prob1, Rodas5());
sol = solve(prob1, Rodas5(autodiff = AutoFiniteDiff()));
