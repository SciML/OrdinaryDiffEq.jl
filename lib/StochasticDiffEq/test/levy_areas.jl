using StochasticDiffEq, DiffEqNoiseProcess, Test, Random, LinearAlgebra
using LevyArea, Statistics

seed = 113213898
Random.seed!(seed)

m = 10
W = WienerProcess(0.0, zeros(m), nothing)

dt = 0.1
calculate_step!(W, dt, nothing, nothing)

for i in 1:10
    accept_step!(W, dt, nothing, nothing)
end

# LevyArea.jl algs: ITER_INT_ALGS

@testset "diagonal noise tests" begin
    true_diag = 1 // 2 .* W.dW .* W.dW

    Jdiag = StochasticDiffEq.JDiagonal_iip(W.dW)
    StochasticDiffEq.get_iterated_I!(dt, W.dW, W.dZ, Jdiag, nothing, 1)
    Jdiagoop = StochasticDiffEq.JDiagonal_oop()

    @test StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, Jdiagoop, nothing, 1) == Jdiag.J
    @test Jdiag.J == true_diag

    for alg in LevyArea.ITER_INT_ALGS
        @test diag(StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, alg, nothing, 1)) ==
            true_diag
    end
end

# Commutative noise tests
@testset "Commutative noise tests" begin
    true_commute = 1 // 2 .* W.dW .* W.dW'

    Jcommute = StochasticDiffEq.JCommute_iip(W.dW)
    StochasticDiffEq.get_iterated_I!(dt, W.dW, W.dZ, Jcommute, nothing, 1)
    Jcommuteoop = StochasticDiffEq.JCommute_oop()

    @test StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, Jcommuteoop, nothing, 1) ==
        Jcommute.J
    @test Jcommute.J == true_commute

    for alg in LevyArea.ITER_INT_ALGS
        @test diag(StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, alg, nothing, 1)) ==
            diag(true_commute)
    end
end

# moment conditions
# E1(I_{j1} I_{j2}) = Δ δ_{j1,j2}
# E2(I_{j1, j2} I_{j1, j2}) = 1/2 Δ^2
# E3(I_{j1} I_{j2} I_{j3, j4}) = {Δ^2 if j1=j2=j3=j4, 1/2 Δ^2 if j3!=j4 and j1=j3, j2=j4 or j1=j4, j2=j3, 0 otherwise}
function test_moments(m, alg, Δ, samples, p = nothing)
    # generate new random dW
    W = WienerProcess!(0.0, zeros(m), nothing)
    E1 = false .* vec(W.dW) .* vec(W.dW)'
    E2 = zero(E1)
    E3 = zero(E1)
    tmp = zero(E1)
    for _ in 1:samples
        calculate_step!(W, Δ, nothing, nothing)
        accept_step!(W, Δ, nothing, nothing)
        mul!(tmp, vec(W.dW), vec(W.dW)')

        I = StochasticDiffEq.get_iterated_I(Δ, W.dW, W.dZ, alg, p, 1)
        @. E1 = E1 + tmp
        @. E2 = E2 + I * I
        @. E3 = E3 + tmp * I
    end
    @. E1 = E1 / samples
    @. E2 = E2 / samples
    @. E3 = E3 / samples

    @test maximum(abs.(E1 - Δ * I)) < 2.0e-1
    @test maximum(abs.(E2 - Δ^2 * (zero(E2) .+ 1) / 2)) < 1.0e-2
    return @test maximum(abs.(E3 - Δ^2 * (one(E2) .+ 1) / 2)) < 1.0e-2
end

"""
Exercise 2.3.1 from
Kloeden, P. E., Platen, E., & Schurz, H. Numerical solution of SDE through computer
experiments. Springer Science & Business Media. (2012)
"""

function test_compare_sample_mean_and_var(alg, Δ, m, samples = Int(1.0e6), p = Int(1.0e2))
    W = WienerProcess(0.0, zeros(m), nothing)
    calculate_step!(W, dt, nothing, nothing)
    for i in 1:10
        accept_step!(W, dt, nothing, nothing)
    end

    xs = [StochasticDiffEq.get_iterated_I(Δ, W.dW, W.dZ, alg, p, 1)[1, 2]]
    coms = [1 // 2 * W.dW[1] * W.dW[2]]

    for i in 1:samples
        calculate_step!(W, Δ, nothing, nothing)
        accept_step!(W, Δ, nothing, nothing)
        com = 1 // 2 * W.dW[1] * W.dW[2]
        x = StochasticDiffEq.get_iterated_I(Δ, W.dW, W.dZ, alg, p, 1)[1, 2]

        push!(xs, x)
        push!(coms, com)
    end

    @test mean(xs) ≈ 0 atol = 1.0e-2
    @test mean(coms) ≈ 0 atol = 1.0e-2
    @test var(xs) ≈ 1 // 2 * Δ^2 rtol = 1.0e-2
    return @test var(coms) ≈ 1 // 4 * Δ^2 rtol = 1.0e-2
end

@testset "General noise tests" begin
    true_commute = 1 // 2 .* W.dW .* W.dW'
    samples = Int(1.0e4)
    for alg in LevyArea.ITER_INT_ALGS
        # Test the relations given in Wiktorsson paper Eq.(2.1)
        Random.seed!(seed)
        I = StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, alg, 1, 1)
        Random.seed!(seed)
        A = LevyArea.levyarea(W.dW / √dt, 1, alg)
        @test dt * A + true_commute == I # because of correction with \mu term
        @test diag(A) == diag(zero(A))
        @test diag(true_commute) == diag(I)
        @test I + I' ≈ 2 * true_commute atol = 1.0e-12
        @test A ≈ -A' atol = 1.0e-12

        # test moment conditions
        test_moments(m, alg, dt, samples)

        # test other StatsJ
        @time test_compare_sample_mean_and_var(alg, 1.0, 2)
    end
end

@testset "Simulate non-com. SDEs tests" begin
    u₀ = [0.0, 0.0, 0.0]
    f(u, p, t) = [0.0, 0.0, 0.0]
    g(u, p, t) = [1.0 0.0; 0.0 1.0; 0.0 u[1]]
    f!(du, u, p, t) = du .*= false
    function g!(du, u, p, t)
        du[1, 1] = 1.0
        du[2, 2] = 1.0
        du[3, 2] = u[1]
        return nothing
    end
    dt = 1 // 2^8
    tspan = (0.0, 1.0)
    prob = SDEProblem(f, g, u₀, (0.0, 1.0); noise_rate_prototype = zeros(3, 2))
    prob! = SDEProblem(f!, g!, u₀, (0.0, 1.0); noise_rate_prototype = zeros(3, 2))

    sol = solve(prob, RKMilGeneral(; ii_approx = IICommutative()), dt = dt, adaptive = false)
    @test sol.retcode == ReturnCode.Success
    sol = solve(prob, RKMilGeneral(; ii_approx = IILevyArea()), dt = dt, adaptive = false)
    @test sol.retcode == ReturnCode.Success
    sol = solve(prob, RKMilGeneral(ii_approx = Fourier()), dt = dt, adaptive = false)
    @test sol.retcode == ReturnCode.Success

    sol = solve(prob!, RKMilGeneral(; ii_approx = IICommutative()), dt = dt, adaptive = false)
    @test sol.retcode == ReturnCode.Success
    sol = solve(prob!, RKMilGeneral(; ii_approx = IILevyArea()), dt = dt, adaptive = false)
    @test sol.retcode == ReturnCode.Success
    sol = solve(prob!, RKMilGeneral(ii_approx = Fourier()), dt = dt, adaptive = false)
    @test sol.retcode == ReturnCode.Success
end
