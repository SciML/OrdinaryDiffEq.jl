using OrdinaryDiffEqNewmark, Test, RecursiveArrayTools, DiffEqDevTools
using SafeTestsets

# Newmark methods with harmonic oscillator
@testset "Harmonic Oscillator" begin
    u0 = fill(0.0, 2)
    v0 = ones(2)
    function f1_harmonic!(dv, v, u, p, t)
        dv .= -u
    end
    function harmonic_jac(J, v, u, weights, p, t)
        J[1, 1] = weights[1] * (0.0) + weights[2] * (-1.0)
        J[1, 2] = weights[1] * (0.0) + weights[2] * (0.0)
        J[2, 2] = weights[1] * (0.0) + weights[2] * (-1.0)
        J[2, 1] = weights[1] * (0.0) + weights[2] * (0.0)
    end
    function f2_harmonic!(du, v, u, p, t)
        du .= v
    end
    function harmonic_analytic(y0, p, x)
        v0, u0 = y0.x
        ArrayPartition(-u0 * sin(x) + v0 * cos(x), u0 * cos(x) + v0 * sin(x))
    end

    ff_harmonic! = DynamicalODEFunction(
        f1_harmonic!, f2_harmonic!; analytic = harmonic_analytic)
    prob = DynamicalODEProblem(ff_harmonic!, v0, u0, (0.0, 5.0))
    dts = 1.0 ./ 2.0 .^ (5:-1:0)

    sim = test_convergence(dts, prob, NewmarkBeta(), dense_errors = true)
    @test sim.𝒪est[:l2]≈2 rtol=1e-1

    function f1_harmonic(v, u, p, t)
        -u
    end
    function f2_harmonic(v, u, p, t)
        v
    end

    ff_harmonic = DynamicalODEFunction(
        f1_harmonic, f2_harmonic; analytic = harmonic_analytic)
    prob = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
    dts = 1.0 ./ 2.0 .^ (5:-1:0)

    sim = test_convergence(dts, prob, NewmarkBeta(), dense_errors = true)
    @test sim.𝒪est[:l2]≈2 rtol=1e-1
end

# Newmark methods with damped oscillator
@testset "Damped Oscillator" begin
    function damped_oscillator!(a, v, u, p, t)
        @. a = -u - 0.5 * v
        return nothing
    end
    function damped_jac(J, v, u, weights, p, t)
        J[1, 1] = weights[1] * (-0.5) + weights[2] * (-1.0)
    end
    function damped_oscillator_analytic(du0_u0, p, t)
        ArrayPartition(
            [
                exp(-t / 4) / 15 * (15 * du0_u0[1] * cos(sqrt(15) * t / 4) -
                 sqrt(15) * (du0_u0[1] + 4 * du0_u0[2]) * sin(sqrt(15) * t / 4))
            ], # du
            [
                exp(-t / 4) / 15 * (15 * du0_u0[2] * cos(sqrt(15) * t / 4) +
                 sqrt(15) * (4 * du0_u0[1] + du0_u0[2]) * sin(sqrt(15) * t / 4))
            ]
        )
    end
    ff_harmonic_damped! = DynamicalODEFunction(
        damped_oscillator!,
        (du, v, u, p, t) -> du = v;
        analytic = damped_oscillator_analytic
    )

    prob = DynamicalODEProblem(ff_harmonic_damped!, [0.0], [1.0], (0.0, 10.0))
    dts = 1.0 ./ 2.0 .^ (5:-1:0)

    sim = test_convergence(dts, prob, NewmarkBeta(), dense_errors = true)
    @test sim.𝒪est[:l2]≈2 rtol=1e-1

    function damped_oscillator(v, u, p, t)
        -u - 0.5 * v
    end
    ff_harmonic_damped = DynamicalODEFunction(
        damped_oscillator,
        (v, u, p, t) -> v;
        analytic = damped_oscillator_analytic
    )

    prob = DynamicalODEProblem(ff_harmonic_damped, [0.0], [1.0], (0.0, 10.0))
    dts = 1.0 ./ 2.0 .^ (5:-1:0)

    sim = test_convergence(dts, prob, NewmarkBeta(), dense_errors = true)
    @test sim.𝒪est[:l2]≈2 rtol=1e-1
end

# Only run JET tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
