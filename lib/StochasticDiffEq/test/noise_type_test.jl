using StochasticDiffEq, SparseArrays, Test, DiffEqNoiseProcess
using SciMLBase: RODEAliasSpecifier

f(du, u, p, t) = (du .= 1.01u)
function g(du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[2]
    du[2, 1] = 1.2u[1]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    return du[2, 4] = 1.8u[2]
end
prob = SDEProblem(f, g, ones(2), (0.0, 1.0), noise_rate_prototype = zeros(2, 4))

sol = solve(prob, EM(), dt = 1 / 1000)

f(du, u, p, t) = (du .= 1.01u)
function g(du, u, p, t)
    du[1, 1] = 0.3
    du[1, 2] = 0.6
    du[1, 3] = 0.9
    du[1, 4] = 0.12
    du[2, 1] = 1.2
    du[2, 2] = 0.2
    du[2, 3] = 0.3
    return du[2, 4] = 1.8
end
prob = SDEProblem(f, g, ones(2), (0.0, 1.0), noise_rate_prototype = zeros(2, 4))

sol = solve(prob, SRA1())
sol = solve(prob, SRA2())
sol = solve(prob, SRA3())
sol = solve(prob, SOSRA())
sol = solve(prob, SOSRA2())
sol = solve(prob, SRA())

@test length(sol.W.u[1]) == 4

f(du, u, p, t) = (du .= 1.01u)
g(du, u, p, t) = (du .= 0.1)
Z = WienerProcess(0.0, [0.0])
prob = SDEProblem(f, g, [1.0], (0.0, 1.0), noise = Z)

sol = solve(prob, EM(), dt = 1 / 100)

@test sol.W == prob.noise
@test objectid(prob.noise) != objectid(sol.W)
@test objectid(prob.noise.u) == objectid(prob.noise.W) != objectid(sol.W.W) ==
    objectid(sol.W.u)

sol = solve(prob, EM(), dt = 1 / 1000, alias = RODEAliasSpecifier(alias_noise = false))

@test sol.W == prob.noise
@test objectid(prob.noise) == objectid(sol.W)
@test objectid(prob.noise.u) == objectid(prob.noise.W) == objectid(sol.W.W) ==
    objectid(sol.W.u)

sol = solve(prob, EM(), dt = 1 / 1000, alias = RODEAliasSpecifier(alias_noise = true))

@test sol.W == prob.noise
@test objectid(prob.noise) != objectid(sol.W)
@test objectid(prob.noise.u) == objectid(prob.noise.W) != objectid(sol.W.W) ==
    objectid(sol.W.u)

function g(du, u, p, t)
    @test du isa SparseMatrixCSC
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[2]
    du[2, 1] = 1.2u[1]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    return du[2, 4] = 1.8u[2]
end
prob = SDEProblem(f, g, ones(2), (0.0, 1.0), noise_rate_prototype = sprand(2, 4, 1.0))

sol = solve(prob, EM(), dt = 1 / 1000)
@test length(sol.W.u[1]) == 4

sol2 = solve(prob, EM(), dt = 1 / 1000)
@test sol.W.curt ≈ sol2.W.curt ≈ 1.0

ff = (u, p, t) -> exp(t)
W = NoiseFunction(0.0, ff)
drift(u, p, t) = u
vol(u, p, t) = u
dt = 0.01
tspan = (0.0, 1.0)
u0 = 0.0
prob = SDEProblem(drift, vol, u0, tspan, noise = W)
sol = solve(prob, EM(), dt = 0.1)
@test sol.W.curt ≈ last(tspan)
sol2 = solve(prob, EM(), dt = 0.1)
@test sol2.W.curt ≈ last(tspan)
tspan = (0.0, 2.0)
prob = SDEProblem(drift, vol, u0, tspan, noise = W)
sol = solve(prob, EM(), dt = 0.01)
@test sol.W.curt ≈ last(tspan)

@test typeof(sol.W) == typeof(prob.noise) && prob.noise isa NoiseFunction
@test objectid(prob.noise) != objectid(sol.W)

sol = solve(prob, EM(), dt = 1 / 1000, alias = RODEAliasSpecifier(alias_noise = false))
@test objectid(prob.noise) == objectid(sol.W)

sol = solve(prob, EM(), dt = 0.01, alias = RODEAliasSpecifier(alias_noise = true))
@test sol.W.curt ≈ last(tspan)

@test typeof(sol.W) == typeof(prob.noise) && prob.noise isa NoiseFunction
@test objectid(prob.noise) != objectid(sol.W)

sol = solve(prob, EM(), dt = 0.01, alias = RODEAliasSpecifier(alias_noise = false))
@test sol.W.curt ≈ last(tspan)

@test typeof(sol.W) == typeof(prob.noise) && prob.noise isa NoiseFunction
@test objectid(prob.noise) == objectid(sol.W)
