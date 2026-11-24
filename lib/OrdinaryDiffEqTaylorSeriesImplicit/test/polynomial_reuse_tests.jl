import OrdinaryDiffEqTaylorSeries

@testset "Implicit Taylor polynomial reuse" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit
    ODTS = OrdinaryDiffEqTaylorSeries

    @testset "Bundle preparation shares a jet without mixing parameters" begin
        calls = Ref(0)
        function counted!(du, u, p, t)
            calls[] += 1
            du .= p .* u .+ t
            return nothing
        end
        low, high = (0.0, 1.0), (0.0, 0.5, -0.2, 0.1)
        specs = (
            (coefficient_sets = (low,), coefficient_jacobians = 0),
            (coefficient_sets = (high, low), coefficient_jacobians = 3),
            (coefficient_sets = (high, low), coefficient_jacobians = 1),
        )
        for p in (-2.0, -3.0)
            before = calls[]
            bundles = ODTS.build_polynomial_bundles(counted!, Val(true), p, specs, 2)
            @test !bundles[1].constant_jacobian
            @test bundles[2].constant_jacobian && bundles[3].constant_jacobian
            @test calls[] == before + 3
            @test bundles === ODTS.build_polynomial_bundles(counted!, Val(true), p, specs, 2)
            @test calls[] == before + 3 # All cache hits: no symbolic preparation.
            for u in ([0.4, 0.2], ComplexF64[0.4 + 0.1im, 0.2 - 0.05im])
                t, h = 0.3, 0.07
                c1 = p .* u .+ t
                c2 = (p .* c1 .+ 1) ./ 2
                c3 = p .* c2 ./ 3
                expected = high[2] * h .* c1 .+ high[3] * h^2 .* c2 .+ high[4] * h^3 .* c3
                a, b = similar(u), similar(u)
                J = fill(eltype(u)(NaN), 2, 2, 3)
                J1 = fill(eltype(u)(NaN), 2, 2, 1)
                bundles[1](a, u, t, h)
                @test a ≈ h .* c1
                bundles[2](a, b, J, u, t, h)
                @test a ≈ expected
                @test b ≈ h .* c1
                for k in 1:3
                    @test J[:, :, k] ≈ (p^k / factorial(k)) .* Matrix{eltype(u)}(I, 2, 2)
                end
                bundles[3](a, b, J1, u, t, h)
                @test a ≈ expected
                @test J1[:, :, 1] == J[:, :, 1]
            end
        end
    end

    @testset "Bundle cache matches parameter values across group and single lookups" begin
        calls = Ref(0)
        function parameterized!(du, u, p, t)
            calls[] += 1
            du .= p[1] .* u .+ p[2] * t
            return nothing
        end
        low, high = (0.0, 1.0), (0.0, 0.5, -0.2, 0.1)
        specs = (
            (coefficient_sets = (low,), coefficient_jacobians = 0),
            (coefficient_sets = (high,), coefficient_jacobians = 0),
        )
        pa, pb = [-2.0, 1.0], [-3.0, 2.0]
        a = ODTS.build_polynomial_bundles(parameterized!, Val(true), pa, specs, 2)
        b = ODTS.build_polynomial_bundles(parameterized!, Val(true), pb, specs, 2)
        @test a !== b
        before = calls[]
        # Preserve value equality, including equal parameters with another
        # container element type; neither identity nor type is a cache key.
        for p in (copy(pa), Float32.(pa))
            @test a === ODTS.build_polynomial_bundles(parameterized!, Val(true), p, specs, 2)
            @test a[1] === ODTS.build_polynomial_bundle(
                parameterized!, Val(true), p, (low,), 2
            )
        end
        @test b === ODTS.build_polynomial_bundles(parameterized!, Val(true), copy(pb), specs, 2)
        @test b[2] === ODTS.build_polynomial_bundle(
            parameterized!, Val(true), copy(pb), (high,), 2
        )
        @test calls[] == before
        u, t, h = [0.4, 0.2], 0.3, 0.07
        for (bundles, p) in ((a, pa), (b, pb))
            out = similar(u)
            c1 = p[1] .* u .+ p[2] * t
            c2 = (p[1] .* c1 .+ p[2]) ./ 2
            c3 = p[1] .* c2 ./ 3
            bundles[1](out, u, t, h)
            @test out ≈ h .* c1
            bundles[2](out, u, t, h)
            @test out ≈ high[2] * h .* c1 .+ high[3] * h^2 .* c2 .+ high[4] * h^3 .* c3
        end
        # The same function and parameters at another state size use another
        # kernel group without displacing the original cached group.
        wider = ODTS.build_polynomial_bundles(parameterized!, Val(true), pa, specs, 3)
        @test wider !== a
        out = zeros(3)
        wider[1](out, [0.4, 0.2, 0.1], t, h)
        @test out ≈ h .* (pa[1] .* [0.4, 0.2, 0.1] .+ pa[2] * t)
        before = calls[]
        @test a === ODTS.build_polynomial_bundles(parameterized!, Val(true), pa, specs, 2)
        @test calls[] == before
    end

    @testset "Joint values and Jacobians" begin
        function rate!(du, u, p, t)
            du[1] = t * u[1] + u[2]^2
            du[2] = p * u[2]
            return nothing
        end
        f = ODEFunction(rate!)
        coefficients = (
            (0.0, 0.5, 1 / 6), (0.0, -0.5, 1 / 6),
            (0.0, 1 / 3), (0.0, -2 / 3, 1 / 3),
        )
        joint = ODTS.build_polynomial_bundle(f, -2.0, coefficients, 2)
        @test joint === ODTS.build_polynomial_bundle(f, -2.0, coefficients, 2)
        pair = ODTS.build_polynomial_bundle(f, -2.0, coefficients[1:2], 2)
        with_jacobian = ODTS.build_polynomial_bundle(
            f, -2.0, (coefficients[2],), 2; jacobian = true
        )
        for u in ([0.4, 0.2], ComplexF64[0.4 + 0.1im, 0.2 - 0.05im])
            t, h = 0.3, 0.07
            outputs = ntuple(_ -> fill!(similar(u), NaN), 4)
            joint(outputs..., u, t, h)
            expected = map(coefficients) do c
                polynomials, _ = ODTS.build_polynomial(f, -2.0, c, 2)
                value = similar(u)
                polynomials[2](value, u, t, h)
                value
            end
            for k in 1:4
                @test outputs[k] ≈ expected[k]
            end
            pair(outputs[1], outputs[2], u, t, h)
            @test outputs[1] ≈ expected[1]
            @test outputs[2] ≈ expected[2]
            J = fill(eltype(u)(NaN), 2, 2)
            with_jacobian(outputs[2], J, u, t, h)
            @test outputs[2] ≈ expected[2]
            @test J[2, 1] == 0 # Overwrite zeros too, including after a prior LU.
            if eltype(u) <: Real
                polynomial = ODTS.build_polynomial(f, -2.0, coefficients[2], 2)[1][1]
                reference = ForwardDiff.jacobian(v -> polynomial(v, t, h), u)
                @test J ≈ reference
            else
                jacobian = ODTS.build_polynomial(f, -2.0, coefficients[2], 2)[2][2]
                reference = similar(J)
                jacobian(reference, u, t, h)
                @test J ≈ reference
            end
        end
    end

    @testset "Step and stage agree" begin
        rate!(du, u, p, t) = (du .= t .* u .+ u .^ 2; nothing)
        problem = ODEProblem(rate!, [0.2, 0.3], (0.0, 0.2))
        for alg in (TaylorGaussA4(), TaylorRadauLA3())
            options = (dt = 0.05, abstol = [1.0e-6, 2.0e-6], reltol = 1.0e-3)
            together, separate = init(problem, alg; options...), init(problem, alg; options...)
            high = similar(problem.u0)
            @test ODTI.implicit_taylor_step!(
                together, together.cache,
                0.0, 0.05, problem.u0, together.u
            )
            @test ODTI.implicit_taylor_stage!(
                separate, separate.cache,
                0.0, 0.05, problem.u0, high
            )
            @test together.u ≈ high rtol = 1.0e-13
            # The embedded estimate is formed from the converged state and the
            # already factorized Newton matrix.
            EEst = OrdinaryDiffEqCore.get_EEst(together)
            @test isfinite(EEst) && EEst > 0
        end
    end

    @testset "Convergence checks use the RHS Jacobian without step-size scaling" begin
        rate!(du, u, p, t) = begin
            du[1] = t * u[1] + u[2]^2
            du[2] = -2u[2]
            nothing
        end
        for u in ([0.4, 0.2], ComplexF64[0.4 + 0.1im, 0.2 - 0.05im])
            problem = ODEProblem(rate!, u, (0.0, 1.0))
            i = init(problem, TaylorGaussA4(); dt = 0.1)
            c = i.cache
            t = 0.3
            expected = [t 2u[2]; zero(eltype(u)) -2]
            for h in (0.0, 0.07, 0.7)
                c.polynomial_check(c.tmp, c.tmp_embedded, c.Jcheck, u, t, h)
                @test c.Mtmp ≈ expected
                first, second = similar(u), similar(u)
                c.polynomial_iter(first, second, u, t, h)
                @test c.tmp ≈ first
                @test c.tmp_embedded ≈ second
            end
        end
    end

    @testset "Matrix reuse remains sensitive to slow rows" begin
        # A large stiff row must not hide a change in a slow row, including
        # rows in the final partial block of the column-wise implementation.
        for T in (Float64, ComplexF64), n in (8, 33, 65)
            kept = zeros(T, n, n)
            for i in 1:n
                kept[i, i] = i == 1 ? 1.0e6 : 1.0
            end
            norms = zeros(n > 32 ? n : 0)
            changes = similar(norms)
            ODTI.matrix_row_norms!(norms, kept)
            current = copy(kept)
            @test !ODTI.rows_changed(current, kept, 0.05, norms, changes)
            current[n, 1] = T <: Complex ? 0.04im : 0.04
            @test !ODTI.rows_changed(current, kept, 0.05, norms, changes)
            current[n, 1] = T <: Complex ? 0.06im : 0.06
            @test ODTI.rows_changed(current, kept, 0.05, norms, changes)
            # The filter compares a negatively scaled Jacobian, not J itself.
            @test ODTI.rows_changed(-current, kept, 0.05, norms, changes, -1)
            copyto!(kept, current)
            ODTI.matrix_row_norms!(norms, kept)
            @test !ODTI.rows_changed(-current, kept, 0.05, norms, changes, -1)
        end
    end

    @testset "Complex center and high orders" begin
        square!(du, u, p, t) = (du .= u .^ 2; nothing)
        problem = ODEProblem(square!, [1.0], (0.0, 0.5))
        for (alg, order) in ((ImplicitTaylor2A4(), 4), (TaylorGaussA6(), 6))
            errors = Float64[]
            for h in (1 / 8, 1 / 16, 1 / 32)
                sol = solve(problem, alg; dt = h, adaptive = false)
                @test SciMLBase.successful_retcode(sol)
                push!(errors, abs(only(sol.u[end]) - 2))
            end
            @test abs(log2(errors[1] / errors[2]) - order) < 0.5
        end
    end
end
