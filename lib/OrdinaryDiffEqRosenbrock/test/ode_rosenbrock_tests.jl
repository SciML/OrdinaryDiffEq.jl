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
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock23())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, Rosenbrock23())
    @test length(sol) < 20
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
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    prob = prob_ode_bigfloat2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock23(linsolve = QRFactorization()))
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, Rosenbrock23(linsolve = QRFactorization()))
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)


    ### ROS34PW3
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW3(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, ROS34PW3())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW3(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, ROS34PW3())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)


    ### ROK4a
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROK4a())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, ROK4a())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROK4a())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, ROK4a())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)


    ### Rodas23W, Rodas3P

    println("Rodas23W")

    prob = prob_ode_linear

    dts = (1 / 2) .^ (6:-1:3)
    sim = test_convergence(dts, prob, Rodas23W(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:L2] ≈ 2 atol = testTol

    sol = solve(prob, Rodas23W())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas23W(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 2 atol = testTol
    @test sim.𝒪est[:L2] ≈ 2 atol = testTol

    sol = solve(prob, Rodas23W())
    @test length(sol) < 20
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
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)
    end


    println("Rodas5P")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5P())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5P())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    println("Rodas5Pe")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5Pe(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5Pe())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5Pe(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5Pe())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    println("Rodas6P")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas6P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 6 atol = testTol

    sol = solve(prob, Rodas6P())
    @test length(sol) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas6P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2] ≈ 6 atol = testTol

    sol = solve(prob, Rodas6P())
    @test length(sol) < 20
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
        @test length(sol) < 20
        @test SciMLBase.successful_retcode(sol)

        prob = prob_ode_2Dlinear

        sim = test_convergence(
            dts, prob,
            Rodas5P(autodiff = AutoEnzyme(mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const)),
            dense_errors = true
        )
        #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
        @test sim.𝒪est[:L2] ≈ 5 atol = testTol

        sim = test_convergence(
            dts, prob,
            Rodas5P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            ),
            dense_errors = true
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
            dense_errors = true
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
        @test length(sol) < 20
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
