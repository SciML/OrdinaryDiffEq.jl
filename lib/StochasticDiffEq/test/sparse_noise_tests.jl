using StochasticDiffEq, SparseArrays, Test, Random

Random.seed!(42)

# 4-state system driven by 2 independent Wiener processes.
# Only states 1 and 3 are coupled to noise; states 2 and 4 are purely deterministic.
# The noise matrix is 4×2 sparse: g[1,1] = 0.1*u[1], g[3,2] = 0.1*u[3], rest zero.

function f_sparse!(du, u, p, t)
    return du .= -u
end

function g_sparse!(du, u, p, t)
    du[1, 1] = 0.1 * u[1]
    return du[3, 2] = 0.1 * u[3]
end

u0 = ones(4)
tspan = (0.0, 1.0)

noise_prototype = spzeros(4, 2)
noise_prototype[1, 1] = 1.0
noise_prototype[3, 2] = 1.0

prob = SDEProblem(f_sparse!, g_sparse!, u0, tspan, noise_rate_prototype = noise_prototype)

@testset "Sparse noise_rate_prototype — EM IIP" begin
    sol = solve(prob, EM(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
    @test length(sol) > 1
    # States 2 and 4 have zero noise rows → evolve purely as du/dt = -u → u(1) = exp(-1)
    @test sol.u[end][2] ≈ exp(-1.0) rtol = 0.02
    @test sol.u[end][4] ≈ exp(-1.0) rtol = 0.02
end

@testset "Sparse noise_rate_prototype — EulerHeun IIP" begin
    sol = solve(prob, EulerHeun(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
    @test length(sol) > 1
    @test sol.u[end][2] ≈ exp(-1.0) rtol = 0.02
    @test sol.u[end][4] ≈ exp(-1.0) rtol = 0.02
end

# OOP variant
function f_sparse_oop(u, p, t)
    return -u
end

function g_sparse_oop(u, p, t)
    du = spzeros(4, 2)
    du[1, 1] = 0.1 * u[1]
    du[3, 2] = 0.1 * u[3]
    return du
end

prob_oop = SDEProblem(f_sparse_oop, g_sparse_oop, u0, tspan, noise_rate_prototype = noise_prototype)

@testset "Sparse noise_rate_prototype — EM OOP" begin
    sol = solve(prob_oop, EM(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end][2] ≈ exp(-1.0) rtol = 0.02
    @test sol.u[end][4] ≈ exp(-1.0) rtol = 0.02
end
