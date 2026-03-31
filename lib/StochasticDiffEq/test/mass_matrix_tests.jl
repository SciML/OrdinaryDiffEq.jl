using StochasticDiffEq, Test, LinearAlgebra, Random

a = ones(3)
t = 0.4
a .- t

const mm_A = [
    -2.0 1 4
    4 -2 1
    2 1 3
]
const mm_b = mm_A * ones(3)
function mm_f(du, u, p, t)
    mul!(du, mm_A, u)
    tmp = t * mm_b
    return @. du += tmp
end
function mm_analytic(u0, p, t, W)
    return 2.0 .* ones(3) .* exp.(t) .- t .- 1
end
function mm_g(du, u, p, t)
    return @. du = u + t
end
function g!(du, u, p, t)
    return @. du = 0.0
end

prob2 = SDEProblem(SDEFunction(mm_g, g!; analytic = mm_analytic), ones(3), (0.0, 1.0))
prob = SDEProblem(
    SDEFunction(mm_f, g!; analytic = mm_analytic, mass_matrix = mm_A),
    ones(3), (0.0, 1.0)
)

sol = solve(prob, ImplicitRKMil(theta = 1), dt = 0.01, adaptive = false)
sol2 = solve(prob2, ImplicitRKMil(theta = 1), dt = 0.01, adaptive = false)

@test norm(sol .- sol2) ≈ 0 atol = 1.0e-11

sol = solve(prob, ImplicitEM(theta = 1), dt = 0.01, adaptive = false)
sol2 = solve(prob2, ImplicitEM(theta = 1), dt = 0.01, adaptive = false)

@test norm(sol .- sol2) ≈ 0 atol = 1.0e-11

sol = solve(prob, ImplicitRKMil(symplectic = true), dt = 0.01, adaptive = false)
sol2 = solve(prob2, ImplicitRKMil(symplectic = true), dt = 0.01, adaptive = false)

@test norm(sol .- sol2) ≈ 0 atol = 1.0e-11

sol = solve(prob, ImplicitEM(symplectic = true), dt = 0.01, adaptive = false)
sol2 = solve(prob2, ImplicitEM(symplectic = true), dt = 0.01, adaptive = false)

@test norm(sol .- sol2) ≈ 0 atol = 1.0e-11

function mm_f2(du, u, p, t)
    return mul!(du, mm_A, u)
end
function no_mm_f2(du, u, p, t)
    return du .= u
end
function no_mm_g2(du, u, p, t)
    return du .= u
end
function mm_g2(du, u, p, t)
    return mul!(du, mm_A, u)
end
prob2 = SDEProblem(no_mm_f2, no_mm_g2, ones(3), (0.0, 1.0))
prob = SDEProblem(SDEFunction(mm_f2, no_mm_g2; mass_matrix = mm_A), ones(3), (0.0, 1.0))

Random.seed!(1)
sol = solve(prob, ImplicitEM(theta = 1), dt = 0.01, adaptive = false)
Random.seed!(1)
sol2 = solve(prob2, ImplicitEM(theta = 1), dt = 0.01, adaptive = false)

@test norm(sol .- sol2) ≈ 0 atol = 1.0e-11
