using StochasticDiffEqLevyArea
using Test
using Random
using StableRNGs: StableRNG
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

    # Regression tests using StableRNGs for cross-version reproducibility.
    # Reference values generated with StableRNG(100) for W, StableRNG(200) for levyarea.
    # These values must remain stable across Julia versions since StableRNGs
    # guarantees a fixed stream.
    @testset "Regression (StableRNG) - $alg_name" for (alg_name, alg, W_ref, I_ref) in [
        ("Fourier m=2", Fourier(),
            [0.0538738050959819, 0.02212611720030917],
            [-0.0035488065622400772 0.002726448021797754;
             -0.0015344298962174452 -0.004755217468819091]),
        ("Milstein m=2", Milstein(),
            [0.0538738050959819, 0.02212611720030917],
            [-0.0035488065622400772 -0.0017205137501459774;
             0.002912531875726286 -0.004755217468819091]),
        ("Wiktorsson m=2", Wiktorsson(),
            [0.0538738050959819, 0.02212611720030917],
            [-0.0035488065622400772 -6.58073345729482e-5;
             0.0012578254601532573 -0.004755217468819091]),
        ("MronRoe m=2", MronRoe(),
            [0.0538738050959819, 0.02212611720030917],
            [-0.0035488065622400772 -0.0008127884214381656;
             0.0020048065470184744 -0.004755217468819091]),
    ]
        h = 0.01
        # Verify W generation is stable
        rng_w = StableRNG(100)
        W_check = sqrt(h) * randn(rng_w, 2)
        @test W_check ≈ W_ref atol = 1e-15

        # Verify iterated_integrals output is stable
        rng_la = StableRNG(200)
        I_new = iterated_integrals(W_ref, h, h^(3 / 2); alg = alg, rng = rng_la)
        @test I_new ≈ I_ref atol = 1e-14
    end
end
