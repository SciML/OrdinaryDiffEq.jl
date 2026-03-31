using StochasticDiffEq, DiffEqNoiseProcess, Random, Test
Random.seed!(100)

f(u, p, t, W) = 1.01u .+ 0.87u .* W
u0 = 1.0
tspan = (0.0, 1.0)
prob = RODEProblem(f, u0, tspan)
sol = solve(prob, RandomEM(), dt = 1 / 100, save_noise = true)
prob = RODEProblem(f, u0, tspan, noise = NoiseWrapper(sol.W))
sol2 = solve(prob, RandomHeun(), dt = 1 / 100)
@test abs(sol.u[end] - sol2.u[end]) < 0.1 * abs(sol.u[end])
sol3 = solve(prob, RandomTamedEM(), dt = 1 / 100)
@test abs(sol.u[end] - sol3.u[end]) < 0.1 * abs(sol.u[end])

f(du, u, p, t, W) = (du .= 1.01u .+ 0.87u .* W)
u0 = ones(4)
prob = RODEProblem(f, u0, tspan)
sol = solve(prob, RandomEM(), dt = 1 / 100, save_noise = true)
prob = RODEProblem(f, u0, tspan, noise = NoiseWrapper(sol.W))
sol2 = solve(prob, RandomHeun(), dt = 1 / 100)
@test sum(abs, sol.u[end] - sol2.u[end]) < 0.1 * sum(abs, sol.u[end])
sol3 = solve(prob, RandomEM(), dt = 1 / 100)
@test sum(abs, sol.u[end] - sol3.u[end]) < 0.1 * sum(abs, sol.u[end])

f(u, p, t, W) = 2u * sin(W)
u0 = 1.0
tspan = (0.0, 1.0)
prob = RODEProblem{false}(f, u0, tspan)
sol = solve(prob, RandomEM(), dt = 1 / 100, save_noise = true)
prob = RODEProblem{false}(f, u0, tspan, noise = NoiseWrapper(sol.W))
sol2 = solve(prob, RandomHeun(), dt = 1 / 100)
@test abs(sol.u[end] - sol2.u[end]) < 0.1 * abs(sol.u[end])
sol3 = solve(prob, RandomTamedEM(), dt = 1 / 100)
@test abs(sol.u[end] - sol3.u[end]) < 0.1 * abs(sol.u[end])

function f(du, u, p, t, W)
    du[1] = 0.2u[1] * sin(W[1] - W[2])
    return du[2] = -0.2u[2] * cos(W[1] + W[2])
end
u0 = [1.0; 1.0]
tspan = (0.0, 1.0)
prob = RODEProblem(f, u0, tspan)
sol = solve(prob, RandomEM(), dt = 1 / 100, save_noise = true)
prob = RODEProblem(f, u0, tspan, noise = NoiseWrapper(sol.W))
sol2 = solve(prob, RandomHeun(), dt = 1 / 100)
@test sum(abs, sol.u[end] - sol2.u[end]) < 0.1 * sum(abs, sol.u[end])
sol3 = solve(prob, RandomTamedEM(), dt = 1 / 100)
@test sum(abs, sol.u[end] - sol3.u[end]) < 0.1 * sum(abs, sol.u[end])

function f(du, u, p, t, W)
    du[1] = -0.2W[3] * u[1] * sin(W[1] - W[2])
    return du[2] = -0.2u[2] * cos(W[1] + W[2])
end
u0 = [1.0; 1.0]
tspan = (0.0, 1.0)
prob = RODEProblem(f, u0, tspan, rand_prototype = zeros(3))
sol = solve(prob, RandomEM(), dt = 1 / 100, save_noise = true)
prob = RODEProblem(f, u0, tspan, noise = NoiseWrapper(sol.W))
sol2 = solve(prob, RandomHeun(), dt = 1 / 100)
@test sum(abs, sol.u[end] - sol2.u[end]) < 0.1 * sum(abs, sol.u[end])
sol3 = solve(prob, RandomTamedEM(), dt = 1 / 100)
@test sum(abs, sol.u[end] - sol3.u[end]) < 0.1 * sum(abs, sol3.u[end])
