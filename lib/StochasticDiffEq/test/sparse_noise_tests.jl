using StochasticDiffEq, SparseArrays, Test, Random

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
    Random.seed!(42)
    sol = solve(prob, EM(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
    @test length(sol) > 1
    # Dimensions 2 and 4 have no noise: pure drift exp(-t)
    @test sol.u[end][2] ≈ exp(-1.0) rtol = 0.02
    @test sol.u[end][4] ≈ exp(-1.0) rtol = 0.02
    # Dimensions 1 and 3 receive noise: path must deviate from pure drift
    @test !isapprox(sol.u[end][1], exp(-1.0), rtol = 0.005)
    @test !isapprox(sol.u[end][3], exp(-1.0), rtol = 0.005)
end

@testset "Sparse noise_rate_prototype — EulerHeun IIP" begin
    Random.seed!(42)
    sol = solve(prob, EulerHeun(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
    @test length(sol) > 1
    @test sol.u[end][2] ≈ exp(-1.0) rtol = 0.02
    @test sol.u[end][4] ≈ exp(-1.0) rtol = 0.02
    @test !isapprox(sol.u[end][1], exp(-1.0), rtol = 0.005)
    @test !isapprox(sol.u[end][3], exp(-1.0), rtol = 0.005)
end

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
    Random.seed!(42)
    sol = solve(prob_oop, EM(), dt = 0.01)
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end][2] ≈ exp(-1.0) rtol = 0.02
    @test sol.u[end][4] ≈ exp(-1.0) rtol = 0.02
    @test !isapprox(sol.u[end][1], exp(-1.0), rtol = 0.005)
    @test !isapprox(sol.u[end][3], exp(-1.0), rtol = 0.005)
end
