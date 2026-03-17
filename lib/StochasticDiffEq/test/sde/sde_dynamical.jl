using StochasticDiffEq, DiffEqNoiseProcess, Test, DiffEqDevTools, Random
Random.seed!(1)

f1_harmonic(v, u, p, t) = -u
f2_harmonic(v, u, p, t) = v
g(u, p, t) = 1
γ = 1

@testset "Scalar u" begin
    u0 = 0
    v0 = 1

    ff_harmonic = DynamicalSDEFunction(f1_harmonic, f2_harmonic, g)
    prob1 = DynamicalSDEProblem(ff_harmonic, v0, u0, (0.0, 5.0))

    dts = (1 / 2) .^ (8:-1:4)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1 = analyticless_test_convergence(
        dts, prob1, BAOAB(gamma = γ), (1 / 2)^10;
        trajectories = Int(2.0e2), use_noise_grid = false
    )
    display(sim1.𝒪est)
    @test abs(sim1.𝒪est[:weak_final] - 1) < 0.3
end

@testset "Vector u" begin
    u0 = zeros(2)
    v0 = ones(2)

    f1_harmonic_iip(dv, v, u, p, t) = dv .= f1_harmonic(v, u, p, t)
    f2_harmonic_iip(du, v, u, p, t) = du .= f2_harmonic(v, u, p, t)
    g_iip(du, u, p, t) = du .= g(u, p, t)

    ff_harmonic = DynamicalSDEFunction(f1_harmonic, f2_harmonic, g)
    prob1 = DynamicalSDEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
    sol1 = solve(prob1, BAOAB(gamma = γ); dt = 1 / 10, save_noise = true)

    prob2 = DynamicalSDEProblem(
        f1_harmonic_iip, f2_harmonic_iip, g_iip, v0,
        u0, (0.0, 5.0); noise = NoiseWrapper(sol1.W)
    )
    sol2 = solve(prob2, BAOAB(gamma = γ); dt = 1 / 10)

    @test sol1[:] ≈ sol2[:]

    dts = (1 / 2) .^ (8:-1:4)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1 = analyticless_test_convergence(
        dts, prob1, BAOAB(gamma = γ), (1 / 2)^10;
        trajectories = Int(1.0e3), use_noise_grid = false
    )
    @test abs(sim1.𝒪est[:weak_final] - 1) < 0.3
end

@testset "Scalar u, scale_noise=false" begin
    u0 = 0
    v0 = 1

    ff_harmonic = DynamicalSDEFunction(f1_harmonic, f2_harmonic, g)
    prob1 = DynamicalSDEProblem(ff_harmonic, v0, u0, (0.0, 5.0))

    dts = (1 / 2) .^ (8:-1:4)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1 = analyticless_test_convergence(
        dts, prob1, BAOAB(gamma = γ, scale_noise = false),
        (1 / 2)^10; trajectories = Int(2.0e2), use_noise_grid = false
    )
    display(sim1.𝒪est)
    @test abs(sim1.𝒪est[:weak_final] - 1) < 0.3
end

@testset "Vector u, scale_noise=false" begin
    u0 = zeros(2)
    v0 = ones(2)

    f1_harmonic_iip(dv, v, u, p, t) = dv .= f1_harmonic(v, u, p, t)
    f2_harmonic_iip(du, v, u, p, t) = du .= f2_harmonic(v, u, p, t)
    g_iip(du, u, p, t) = du .= g(u, p, t)

    ff_harmonic = DynamicalSDEFunction(f1_harmonic, f2_harmonic, g)
    prob1 = DynamicalSDEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
    sol1 = solve(prob1, BAOAB(gamma = γ, scale_noise = false); dt = 1 / 10, save_noise = true)

    prob2 = DynamicalSDEProblem(
        f1_harmonic_iip, f2_harmonic_iip, g_iip, v0,
        u0, (0.0, 5.0); noise = NoiseWrapper(sol1.W)
    )
    sol2 = solve(prob2, BAOAB(gamma = γ, scale_noise = false); dt = 1 / 10)

    @test sol1[:] ≈ sol2[:]

    dts = (1 / 2) .^ (8:-1:4)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1 = analyticless_test_convergence(
        dts, prob1, BAOAB(gamma = γ, scale_noise = false),
        (1 / 2)^10; trajectories = Int(1.0e3), use_noise_grid = false
    )
    @test abs(sim1.𝒪est[:weak_final] - 1) < 0.3
end
