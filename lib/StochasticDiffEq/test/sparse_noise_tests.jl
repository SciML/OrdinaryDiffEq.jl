using StochasticDiffEq, SparseArrays, Test, Random

# Four states driven by two independent Wiener processes. g is structurally sparse:
# channel 1 drives state 1, channel 2 drives state 3, states 2 and 4 are pure drift.

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
# The same problem with a dense g prototype: the reference the sparse run has to match.
prob_dense = SDEProblem(
    f_sparse!, g_sparse!, u0, tspan, noise_rate_prototype = Matrix(noise_prototype)
)

@testset "Sparse noise_rate_prototype — dW is dense" begin
    integrator = init(prob, EM(), dt = 0.01)
    @test !issparse(integrator.W.dW)
    @test integrator.W.dW isa Vector{Float64}
end

@testset "Sparse noise_rate_prototype — caches with a dense dZ" begin
    # W2Ito1 sizes dZ densely no matter what dW is and stores both under one cache
    # field type, so a sparse dW leaves the cache constructor unresolvable.
    @test solve(prob, W2Ito1(), dt = 0.01, adaptive = false).retcode == ReturnCode.Success
end

@testset "Sparse noise_rate_prototype — same path as a dense prototype" begin
    # RKMilGeneral draws a long dZ for the Lévy area terms. Sparse and dense arrays of
    # that length take different randn! methods, so the same seed gives a different
    # realisation unless the prototype is dense.
    Random.seed!(42)
    sol_sparse = solve(prob, RKMilGeneral(p = 10), dt = 0.01, adaptive = false)
    Random.seed!(42)
    sol_dense = solve(prob_dense, RKMilGeneral(p = 10), dt = 0.01, adaptive = false)
    @test sol_sparse.u[end] == sol_dense.u[end]
end

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
