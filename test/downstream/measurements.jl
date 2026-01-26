using Measurements, OrdinaryDiffEq, StaticArrays

function lorenz_rule(u, p, t)
    σ, ρ, β = p
    x, y, z = u
    dx = σ * (y - x)
    dy = x * (ρ - z) - y
    dz = x * y - β * z
    return SVector(dx, dy, dz)
end

u₀ = SVector(10.0, 10.0, 10.0)
p₀ = [10, 28, 8 / 3]

prob = ODEProblem(lorenz_rule, u₀, (0.0, 100.0), p₀)
alg = Vern9() # 9-th order adaptive solver
sol = solve(prob; alg = alg)

uunc = SVector(10.0 ± 0.1, 10.0 ± 0.1, 10.0 ± 0.1)
probunc = ODEProblem(lorenz_rule, uunc, (0.0, 10.0), p₀)
solunc = solve(probunc; alg = alg)
