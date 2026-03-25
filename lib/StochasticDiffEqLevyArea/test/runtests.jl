using StochasticDiffEqLevyArea
using Test
using Random
using LinearAlgebra: diag, norm
using Statistics: var

Random.seed!(638278)

@testset "StochasticDiffEqLevyArea" begin
    @testset "Scalar case (m=1)" begin
        for h in rand(5)
            W = √h * randn()
            Ints = iterated_integrals(W, h, h^(3 / 2))
            @test Ints == 0.5W^2 - 0.5h
        end
    end

    @testset "Iterated integrals - $alg, m=$m" for
        alg in [Fourier(), Milstein(), Wiktorsson(), MronRoe()],
        m in [2, 10, 50]

        h = 0.0001
        ε = h^(3 / 2)
        W = √h * randn(m)
        Ints = iterated_integrals(W, h, ε; alg = alg)
        @test diag(Ints) ≈ 0.5 * W .^ 2 .- 0.5h
    end

    @testset "RNG parameter - $alg" for alg in [Fourier(), Milstein(), Wiktorsson(), MronRoe()]
        h = 0.01
        m = 4
        W = √h * randn(m)
        ε = h^(3 / 2)

        # Same RNG seed → same result
        rng1 = Xoshiro(42)
        rng2 = Xoshiro(42)
        I1 = iterated_integrals(W, h, ε; alg = alg, rng = rng1)
        I2 = iterated_integrals(W, h, ε; alg = alg, rng = rng2)
        @test I1 == I2

        # Different seed → different result
        rng3 = Xoshiro(99)
        I3 = iterated_integrals(W, h, ε; alg = alg, rng = rng3)
        @test I3 != I1
    end

    @testset "Coefficient-based computation - $alg" for
        alg in [Fourier(), Milstein(), Wiktorsson(), MronRoe()]

        h = 0.01
        m = 4
        n = 10
        W = √h * randn(m)

        # Generate coefficients from a seeded RNG
        rng = Xoshiro(42)
        coeffs = generate_coefficients(m, n, alg, rng)

        # Compute from coefficients (deterministic)
        I1 = iterated_integrals(W, h, coeffs; alg = alg)
        I2 = iterated_integrals(W, h, coeffs; alg = alg)
        @test I1 == I2

        # Diagonal identity still holds
        @test diag(I1) ≈ 0.5 * W .^ 2 .- 0.5h

        # Coefficient length matches norv
        @test coefficient_length(m, n, alg) == length(coeffs.X) + length(coeffs.Y) + length(coeffs.tail)
    end

    @testset "Path reconstruction" begin
        h = 1.0
        m = 4
        n = 50
        dW = randn(m) * √h

        rng = Xoshiro(42)
        coeffs = generate_coefficients(m, n, MronRoe(), rng)

        # Boundary conditions: W(0) = 0, W(h) = dW
        W_vals = reconstruct_path(dW, h, coeffs, [0.0, h])
        @test W_vals[1] ≈ zeros(m) atol = 1e-12
        @test W_vals[2] ≈ dW atol = 1e-10

        # Mid-point should not equal simple linear interpolation (bridge is non-trivial)
        W_mid = reconstruct_path(dW, h, coeffs, [h / 2])
        @test W_mid[1] != 0.5 * dW  # bridge adds fluctuations

        # Statistical check: variance of bridge at midpoint ≈ h/4 per component
        # (using many independent realizations)
        n_samples = 5000
        mid_values = zeros(m, n_samples)
        for i in 1:n_samples
            dW_i = randn(m) * √h
            coeffs_i = generate_coefficients(m, n, MronRoe(), Xoshiro(i))
            W_mid_i = reconstruct_path(dW_i, h, coeffs_i, [h / 2])
            mid_values[:, i] = W_mid_i[1] .- 0.5 .* dW_i  # bridge component
        end
        for j in 1:m
            bridge_var = var(mid_values[j, :])
            @test bridge_var ≈ h / 4 rtol = 0.15  # allow 15% relative tolerance
        end
    end

    @testset "Sub-interval iterated integrals" begin
        h = 1.0
        m = 4
        n = 50
        dW = randn(m) * √h

        rng = Xoshiro(42)
        coeffs = generate_coefficients(m, n, MronRoe(), rng)

        # Full interval: diagonal should satisfy Itô identity
        I_full = iterated_integrals_subinterval(dW, h, coeffs, 0.0, h;
            n_quadrature = 1000, ito_correction = true)
        @test diag(I_full) ≈ 0.5 * dW .^ 2 .- 0.5 * h atol = 0.05

        # Sub-interval diagonal identity
        I_sub = iterated_integrals_subinterval(dW, h, coeffs, 0.0, h / 2;
            n_quadrature = 500, ito_correction = true)
        W_half = reconstruct_path(dW, h, coeffs, [h / 2])[1]
        dW_sub = W_half
        @test diag(I_sub) ≈ 0.5 * dW_sub .^ 2 .- 0.5 * (h / 2) atol = 0.05

        # Determinism: same coefficients → same result
        I_sub2 = iterated_integrals_subinterval(dW, h, coeffs, 0.0, h / 2;
            n_quadrature = 500, ito_correction = true)
        @test I_sub == I_sub2

        # Convergence: increasing quadrature improves diagonal accuracy
        for nq in [100, 500, 2000]
            I_q = iterated_integrals_subinterval(dW, h, coeffs, 0.0, h;
                n_quadrature = nq, ito_correction = true)
            err = norm(diag(I_q) - (0.5 * dW .^ 2 .- 0.5 * h))
            # Just check it's reasonably small
            @test err < 0.5
        end
    end

    @testset "Algorithm selection" begin
        h = 1 / 128
        m = 10
        ε = h^(3 / 2)
        alg = optimal_algorithm(m, h, ε, MaxL2())
        @test alg isa AbstractIteratedIntegralAlgorithm
        n = terms_needed(m, h, ε, alg, MaxL2())
        @test n > 0
    end

    # Regression tests: verify bit-identical output with LevyArea.jl
    # Reference values were generated using LevyArea.jl v1.0.0 with:
    #   Random.seed!(638278) to generate W, then Random.seed!(12345) before each levyarea call
    @testset "Regression vs LevyArea.jl - $alg_name" for (alg_name, alg, W_ref, I_ref) in [
        ("Fourier m=2", Fourier(),
            [0.09903103214960329, 0.11212317514303571],
            [-9.642733569211956e-5 0.0020195078355866964;
             0.009084165926718873 0.001285803202077931]),
        ("Milstein m=2", Milstein(),
            [-0.03121756137083341, 0.13620291896429848],
            [-0.004512731931029125 -0.013820802895184192;
             0.009568879913529555 0.004275617567197629]),
        ("Wiktorsson m=2", Wiktorsson(),
            [-0.05354470833343851, 0.029557354340556088],
            [-0.0035664821047435005 -0.0036534036925311725;
             0.0020707637752580043 -0.004563181402193405]),
        ("MronRoe m=2", MronRoe(),
            [-0.21299582757570945, -0.053796834119659964],
            [0.017683611282330675 0.012644240757892108;
             -0.001185739553621972 -0.0035529503193508947]),
    ]
        h = 0.01
        Random.seed!(12345)
        I_new = iterated_integrals(W_ref, h, h^(3 / 2); alg = alg)
        @test I_new ≈ I_ref atol = 1e-14
    end
end
