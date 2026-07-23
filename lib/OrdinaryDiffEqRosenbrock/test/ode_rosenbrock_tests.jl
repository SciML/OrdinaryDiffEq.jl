using OrdinaryDiffEqRosenbrock, DiffEqDevTools, Test, LinearAlgebra, LinearSolve, ADTypes
import ODEProblemLibrary: prob_ode_linear,
    prob_ode_2Dlinear,
    prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear
import LinearSolve

if isempty(VERSION.prerelease)
    using Enzyme
end

@testset "Rosenbrock Tests" begin
    ## Breakout these since no other test of their adaptivity

    dts = (1 / 2) .^ (6:-1:3)
    testTol = 0.2

    ### Rosenbrock23()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rosenbrock23())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, Rosenbrock23())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock23())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, Rosenbrock23())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts, prob, Rosenbrock23(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                )
            )
        )
        @test sim.𝒪est[:final] ≈ 2 atol = testTol

        sol = solve(
            prob, Rosenbrock23(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                )
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    prob = prob_ode_bigfloat2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock23(linsolve = QRFactorization()))
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, Rosenbrock23(linsolve = QRFactorization()))
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)


    ### ROS34PW3
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW3(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, ROS34PW3())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW3(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, ROS34PW3())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)


    ### ROK4a
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROK4a())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, ROK4a())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROK4a())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, ROK4a())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)


    ### Rodas23W, Rodas3P

    println("Rodas23W")

    prob = prob_ode_linear

    dts = (1 / 2) .^ (6:-1:3)
    sim = test_convergence(dts, prob, Rodas23W(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:L2] ≈ 2 atol = testTol

    sol = solve(prob, Rodas23W())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas23W(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:L2] ≈ 2 atol = testTol

    sol = solve(prob, Rodas23W())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            Rodas23W(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 2 atol = testTol

        sol = solve(
            prob,
            Rodas23W(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end


    println("Rodas3d")

    # Rodas3d's damping parameter γ = 0.57281606 is a root of the 4th-order
    # linear order condition (an endpoint of the L-stability interval in
    # Hairer & Wanner Table 6.4), so the method superconverges at order 4 on
    # linear problems; the design order 3 is checked on a nonlinear problem below.
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas3d())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Rodas3d())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas3d())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Rodas3d())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = ODEProblem(
        (u, p, t) -> [-2u[1] + u[2]^2, -u[2] + sin(u[1]) + 0.1t],
        [1.0, 0.5], (0.0, 1.0)
    )
    test_setup = Dict(:alg => Rodas4(), :reltol => 1.0e-13, :abstol => 1.0e-13)
    sim = analyticless_test_convergence((1 / 2) .^ (7:-1:4), prob, Rodas3d(), test_setup)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    println("Rodas5P")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    println("Rodas5Pe")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5Pe(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5Pe())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5Pe(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5Pe())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    println("Rodas6P")

    prob = prob_ode_linear

    rodas6p_dts = (1 / 2) .^ (5:-1:2)
    sim = test_convergence(rodas6p_dts, prob, Rodas6P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 6 atol = testTol

    sol = solve(prob, Rodas6P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(rodas6p_dts, prob, Rodas6P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 6 atol = testTol

    sol = solve(prob, Rodas6P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)


    println("Rodas5P Enzyme Forward")

    prob = prob_ode_linear

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts, prob,
            Rodas5P(autodiff = AutoEnzyme(mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const)),
            dense_errors = true
        )
        #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
        @test sim.𝒪est[:L2] ≈ 5 atol = testTol

        sol = solve(
            prob,
            Rodas5P(autodiff = AutoEnzyme(mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const))
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)

        prob = prob_ode_2Dlinear

        sim = test_convergence(
            dts, prob,
            Rodas5P(autodiff = AutoEnzyme(mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const)),
            dense_errors = true
        )
        #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
        @test sim.𝒪est[:L2] ≈ 5 atol = testTol

        # Krylov linear solvers need tight tolerance for the outer-method
        # convergence order to surface; LinearSolve now respects reltol/abstol.
        sim = test_convergence(
            dts, prob,
            Rodas5P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            ),
            dense_errors = true,
            reltol = 1.0e-14, abstol = 1.0e-14
        )
        #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
        @test sim.𝒪est[:L2] ≈ 5 atol = testTol

        sim = test_convergence(
            dts, prob,
            Rodas5P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL_GMRES()
            ),
            dense_errors = true,
            reltol = 1.0e-14, abstol = 1.0e-14
        )
        #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
        @test sim.𝒪est[:L2] ≈ 5 atol = testTol

        sol = solve(
            prob,
            Rodas5P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward),
                    function_annotation = Enzyme.Const
                )
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)


        prob = ODEProblem((u, p, t) -> 0.9u, 0.1, (0.0, 1.0))
        @test_nowarn solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
        @test_nowarn solve(
            prob,
            Rosenbrock23(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward),
                    function_annotation = Enzyme.Const
                )
            )
        )
    end
end

@testset "Convergence with time-dependent matrix-free Jacobian" begin
    time_derivative(du, u, p, t) = (du[1] = t * u[1])
    time_derivative_analytic(u0, p, t) = u0 * exp(t^2 / 2)
    ff_time_derivative = ODEFunction(time_derivative, analytic = time_derivative_analytic)
    prob = ODEProblem(ff_time_derivative, [1.0], (0.0, 1.0))

    dts = (1 / 2) .^ (6:-1:3)
    testTol = 0.2
    # Check convergence of Rodas3 with time-dependent matrix-free Jacobian.
    # Primarily to check that the Jacobian is being updated correctly as t changes.
    sim = test_convergence(
        dts, prob, Rodas3(linsolve = LinearSolve.KrylovJL()),
        reltol = 1.0e-14, abstol = 1.0e-14
    )
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
end

@testset "ADTypes" begin
    for T in [
            Rosenbrock23,
            Rodas23W,
            Rodas5P,
            Rodas5Pe,
            Rodas5Pr,
            Rodas6P,
            ROS34PW3,
            ROK4a,
        ]
        RosenbrockAlgorithm = if T <:
            OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAlgorithm
            OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAlgorithm
        else
            OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
        end

        # Test with AutoForwardDiff
        ad = AutoForwardDiff(; chunksize = 3)
        alg = T(; autodiff = ad)
        @test alg isa RosenbrockAlgorithm
        @test alg.autodiff === ad
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.alg_autodiff(alg) === ad
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.get_chunksize(alg) === Val{3}()

        # chunk_size keyword was removed in v7, test that it errors
        @test_throws MethodError T(; autodiff = ad, chunk_size = Val{4}())

        # Test with AutoFiniteDiff
        ad = AutoFiniteDiff(; fdtype = Val{:central}())
        alg = T(; autodiff = ad)
        @test alg isa RosenbrockAlgorithm
        @test alg.autodiff isa AutoFiniteDiff{Val{:central}}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.alg_autodiff(alg) isa AutoFiniteDiff{Val{:central}}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.get_chunksize(alg) === Val{0}()

        # diff_type keyword was removed in v7, test that it errors
        @test_throws MethodError T(; autodiff = ad, diff_type = Val{:complex}())

        # issue #2613
        f(u, _, _) = -u
        prob = ODEProblem(f, [1.0, 0.0], (0.0, 1.0))
        alg = T(; autodiff = AutoForwardDiff(; chunksize = 1))
        sol = if alg isa OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
            @inferred(solve(prob, alg))
        else
            @inferred(solve(prob, alg; dt = 0.1))
        end
        alg = T(; autodiff = AutoFiniteDiff(; fdtype = Val(:central)))
        sol = if alg isa OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
            @inferred(solve(prob, alg))
        else
            @inferred(solve(prob, alg; dt = 0.1))
        end
    end
end

@testset "Issue #3631: DAE dense output for Rosenbrock methods with empty H" begin
    # Methods without a stiff-aware dense output (empty H matrix, e.g. ROS34PW2,
    # ROS34PW3, Rodas3) stored f(uprev) and f(u) and reused the standard Hermite
    # cubic. For DAE problems those values are residuals on algebraic variables,
    # not derivatives, so Hermite was producing wildly wrong off-knot values.
    # See https://github.com/SciML/OrdinaryDiffEq.jl/issues/3631.
    function f!(du, u, p, t)
        du[1] = 20 * ((t - 1)^(p - 1) + t * (p - 1) * (t - 1)^(p - 2))
        du[2] = u[1] - u[2]
        return nothing
    end
    anasol(t, p) = 10 .* t .* (t .- 1) .^ (p - 1)

    p = 3
    u0 = [0.0, 0.0]
    tspan = (0.0, 1.5)
    M = zeros(2, 2)
    M[1, 1] = 1
    M[1, 2] = 1
    prob = ODEProblem(ODEFunction(f!, mass_matrix = M), u0, tspan, p)

    tt = collect(0:0.01:tspan[2])

    # Methods with their own dense output should match the analytical solution
    # almost exactly; methods that fall back to the generic Hermite path must
    # at least be bounded by O(dt) — Hermite-on-residuals previously produced
    # interpolation errors orders of magnitude larger than the knot errors.
    for (method, expected_interp) in (
            (Rodas4P, 1.0e-10), (ROS34PW2, 5.0e-2),
            (ROS34PW3, 5.0e-2), (Rodas3, 5.0e-2), (Rodas3d, 5.0e-2),
        )
        sol = solve(prob, method(); dense = true, reltol = 1.0e-4, abstol = 1.0e-4)
        err_knot = maximum(abs.(sol[1, :] .- anasol(sol.t, p)))
        err_interp = maximum(abs.(sol(tt; idxs = 1) .- anasol(tt, p)))
        @test err_knot < 1.0e-6
        @test err_interp < expected_interp
        # interp error should not be wildly larger than knot error
        @test err_interp < 100 * max(err_knot, 1.0e-4)
    end
end

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/3721
# OOP static-array mass-matrix solves were failing in two ways:
# 1. With a Float32 state and a Float64 time type (e.g. `dt = 0.1` promoting
#    tspan), calc_W produces a Float64-eltype StaticWOperator that cannot be
#    assigned into JacReuseState's Float32-seeded `cached_W` slot.
# 2. On rejected-step retries the OOP reuse branch rebuilt W via
#    `default_factorize`, producing a `StaticArrays.LU` instead of the
#    `StaticWOperator` type of the `cached_W` slot.
using StaticArrays
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit

@testset "Static-array mass-matrix OOP JacReuseState (#3721)" begin
    function rober_static(u, p, t)
        y₁, y₂, y₃ = u
        k₁, k₂, k₃ = p
        return @SVector [
            -k₁ * y₁ + k₃ * y₂ * y₃,
            k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃,
            y₁ + y₂ + y₃ - 1,
        ]
    end

    # Float64 reference
    M64 = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
    prob64 = ODEProblem(
        ODEFunction(rober_static, mass_matrix = M64),
        @SVector([1.0, 0.0, 0.0]), (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4)
    )
    ref = solve(
        prob64, Rosenbrock23(), abstol = 1.0e-10, reltol = 1.0e-10,
        initializealg = BrownFullBasicInit()
    )

    # Failure mode 1: Float32 state with Float64 dt (promotes the time type).
    M32 = @SMatrix [1.0f0 0.0f0 0.0f0; 0.0f0 1.0f0 0.0f0; 0.0f0 0.0f0 0.0f0]
    prob32 = ODEProblem(
        ODEFunction(rober_static, mass_matrix = M32),
        @SVector([1.0f0, 0.0f0, 0.0f0]), (0.0f0, 1.0f5), (0.04f0, 3.0f7, 1.0f4)
    )
    sol32 = solve(
        prob32, Rosenbrock23(), dt = 0.1, abstol = 1.0f-5, reltol = 1.0f-5,
        initializealg = BrownFullBasicInit()
    )
    @test SciMLBase.successful_retcode(sol32)
    @test norm(sol32.u[end] - ref.u[end]) < 1.0e-4

    # Failure mode 2: homogeneous Float64 at a tolerance that produces step
    # rejections, whose retries take the cached-J reuse branch. nreject > 0 is
    # asserted so this keeps covering that branch if step control changes.
    sol64 = solve(
        prob64, Rosenbrock23(), abstol = 1.0e-5, reltol = 1.0e-5,
        initializealg = BrownFullBasicInit()
    )
    @test SciMLBase.successful_retcode(sol64)
    @test sol64.stats.nreject > 0
    @test norm(sol64.u[end] - ref.u[end]) < 1.0e-4

    # The Rodas-family OOP constant cache goes through the same cached_W
    # seeding path; Rodas23W and ROS34PW2 are W-methods so they additionally
    # exercise the Jacobian-reuse decision logic on the static path.
    for alg in (Rodas4(), Rodas23W(), ROS34PW2())
        sol_rodas = solve(
            prob32, alg, dt = 0.1, abstol = 1.0f-5, reltol = 1.0f-5,
            initializealg = BrownFullBasicInit()
        )
        @test SciMLBase.successful_retcode(sol_rodas)
        @test norm(sol_rodas.u[end] - ref.u[end]) < 1.0e-4
    end
end
