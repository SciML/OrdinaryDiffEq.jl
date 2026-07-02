using Pkg
using OrdinaryDiffEqNewmark, Test, RecursiveArrayTools, DiffEqDevTools
using LinearAlgebra: norm
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    # Newmark methods with harmonic oscillator
    @testset "Newmark Harmonic Oscillator" begin
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
            f1_harmonic!, f2_harmonic!; analytic = harmonic_analytic
        )
        prob = DynamicalODEProblem(ff_harmonic!, v0, u0, (0.0, 5.0))
        dts = 1.0 ./ 2.0 .^ (5:-1:0)

        sim = test_convergence(dts, prob, NewmarkBeta(), dense_errors = true)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1

        function f1_harmonic(v, u, p, t)
            -u
        end
        function f2_harmonic(v, u, p, t)
            v
        end

        ff_harmonic = DynamicalODEFunction(
            f1_harmonic, f2_harmonic; analytic = harmonic_analytic
        )
        prob = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
        dts = 1.0 ./ 2.0 .^ (5:-1:0)

        sim = test_convergence(dts, prob, NewmarkBeta(), dense_errors = true)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
    end

    # Newmark methods with damped oscillator
    @testset "Newmark Damped Oscillator" begin
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
                    exp(-t / 4) / 15 * (
                        15 * du0_u0[1] * cos(sqrt(15) * t / 4) -
                            sqrt(15) * (du0_u0[1] + 4 * du0_u0[2]) * sin(sqrt(15) * t / 4)
                    ),
                ], # du
                [
                    exp(-t / 4) / 15 * (
                        15 * du0_u0[2] * cos(sqrt(15) * t / 4) +
                            sqrt(15) * (4 * du0_u0[1] + du0_u0[2]) * sin(sqrt(15) * t / 4)
                    ),
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
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1

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
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
    end

    @testset "GeneralizedAlpha Harmonic Oscillator" begin
        # Convergence uses saved-step errors (`dense_errors=false`): dense Hermite output is not tailored to GA.
        u0 = fill(0.0, 2)
        v0 = ones(2)
        function f1_harmonic!(dv, v, u, p, t)
            dv .= -u
        end
        function f2_harmonic!(du, v, u, p, t)
            du .= v
        end
        function harmonic_analytic(y0, p, x)
            v0, u0 = y0.x
            ArrayPartition(-u0 * sin(x) + v0 * cos(x), u0 * cos(x) + v0 * sin(x))
        end
        ff_harmonic! = DynamicalODEFunction(f1_harmonic!, f2_harmonic!; analytic = harmonic_analytic)
        prob = DynamicalODEProblem(ff_harmonic!, v0, u0, (0.0, 5.0))
        dts = 1.0 ./ 2.0 .^ (5:-1:0)

        # rho_inf = 1 is equivalent to undamped Newmark, should be 2nd order
        sim = test_convergence(dts, prob, GeneralizedAlpha(rho_inf = 1.0), dense_errors = false)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1

        # HHT constructor, when 0 equal to Newmark
        sim = test_convergence(dts, prob, GeneralizedAlpha(alpha_hht = 0.0), dense_errors = false)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1

        # Direct 4-parameter constructor, same params as rho_inf = 1.0
        sim = test_convergence(dts, prob, GeneralizedAlpha(0.0, 0.0, 0.25, 0.5), dense_errors = false)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1

        # Out-of-place
        function f1_harmonic(v, u, p, t)
            -u
        end
        function f2_harmonic(v, u, p, t)
            v
        end
        ff_harmonic = DynamicalODEFunction(f1_harmonic, f2_harmonic; analytic = harmonic_analytic)
        prob_oop = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))
        sim = test_convergence(dts, prob_oop, GeneralizedAlpha(rho_inf = 1.0), dense_errors = false)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
    end

    @testset "GeneralizedAlpha Damped Oscillator" begin
        function damped_oscillator!(a, v, u, p, t)
            @. a = -u - 0.5 * v
        end
        function damped_oscillator_analytic(du0_u0, p, t)
            ArrayPartition(
                [
                    exp(-t / 4) / 15 * (
                        15 * du0_u0[1] * cos(sqrt(15) * t / 4) -
                            sqrt(15) * (du0_u0[1] + 4 * du0_u0[2]) * sin(sqrt(15) * t / 4)
                    ),
                ],
                [
                    exp(-t / 4) / 15 * (
                        15 * du0_u0[2] * cos(sqrt(15) * t / 4) +
                            sqrt(15) * (4 * du0_u0[1] + du0_u0[2]) * sin(sqrt(15) * t / 4)
                    ),
                ]
            )
        end
        ff_damped! = DynamicalODEFunction(
            damped_oscillator!, (du, v, u, p, t) -> du .= v;
            analytic = damped_oscillator_analytic
        )
        prob = DynamicalODEProblem(ff_damped!, [0.0], [1.0], (0.0, 10.0))
        dts = 1.0 ./ 2.0 .^ (5:-1:0)

        sim = test_convergence(dts, prob, GeneralizedAlpha(rho_inf = 1.0), dense_errors = false)
        @test sim.𝒪est[:l2] ≈ 2 rtol = 1.0e-1
    end

    @testset "GeneralizedAlpha Algorithmic Damping" begin
        # High frequency problem — solution is a fast oscillation on top of a slow one.
        # rho_inf = 0 should damp the high frequency content, rho_inf = 1 should not.
        function f1_highfreq!(dv, v, u, p, t)
            @. dv = -100.0 * u  # fast mode: ω = 10
        end
        function f2_highfreq!(du, v, u, p, t)
            du .= v
        end
        ff_highfreq! = DynamicalODEFunction(f1_highfreq!, f2_highfreq!)
        # initial condition excites only the fast mode
        prob = DynamicalODEProblem(ff_highfreq!, ones(1), zeros(1), (0.0, 1.0))

        sol_damped = solve(prob, GeneralizedAlpha(rho_inf = 0.5), dt = 0.1, adaptive = false)
        sol_undamped = solve(prob, GeneralizedAlpha(rho_inf = 1.0), dt = 0.1, adaptive = false)

        # damped solution should have smaller amplitude at end of integration
        @test norm(sol_damped.u[end]) < norm(sol_undamped.u[end])
    end

    @testset "GeneralizedAlpha Constructor Validation" begin
        # rho_inf out of range
        @test_throws AssertionError GeneralizedAlpha(rho_inf = 1.1)
        @test_throws AssertionError GeneralizedAlpha(rho_inf = -0.1)

        # HHT out of range
        @test_throws AssertionError GeneralizedAlpha(alpha_hht = 0.1)
        @test_throws AssertionError GeneralizedAlpha(alpha_hht = -0.4)

        # direct params violating stability
        @test_throws AssertionError GeneralizedAlpha(0.3, 0.1, 0.25, 0.5)
        @test_throws AssertionError GeneralizedAlpha(0.0, 0.6, 0.25, 0.5)
    end
end

# Run QA tests (AllocCheck, JET) - skip on pre-release Julia
# Allocation tests must run before JET because JET's static analysis
# invalidates compiled code and causes spurious runtime allocations.
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
end
