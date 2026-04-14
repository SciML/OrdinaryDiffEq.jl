using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqRosenbrockTableaus
using DiffEqDevTools, Test, LinearAlgebra, LinearSolve, ADTypes
import ODEProblemLibrary: prob_ode_linear,
    prob_ode_2Dlinear,
    prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear
import LinearSolve

if isempty(VERSION.prerelease)
    using Enzyme
end

@testset "Non-essential Rosenbrock Tests" begin
    dts = (1 / 2) .^ (6:-1:3)
    testTol = 0.2

    ### Rosenbrock32()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rosenbrock32())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, Rosenbrock32())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock32())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, Rosenbrock32())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            Rosenbrock32(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                )
            )
        )
        @test sim.𝒪est[:final] ≈ 3 atol = testTol

        sol = solve(
            prob,
            Rosenbrock32(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                )
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)

        sim = test_convergence(
            dts,
            prob,
            Rosenbrock32(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ), linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 3 atol = testTol

        sol = solve(
            prob,
            Rosenbrock32(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ), linsolve = LinearSolve.KrylovJL()
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end
    ### ROS3P()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3P(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS3P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3P(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS3P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            ROS3P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 3 atol = testTol

        sol = solve(
            prob,
            ROS3P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    ### Rodas3()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas3(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, Rodas3())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas3(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, Rodas3())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            Rodas3(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 3 atol = testTol

        sol = solve(
            prob,
            Rodas3(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    ### ROS2
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, ROS2())
    @test length(sol.t) < 61
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS2())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, ROS2PR())
    @test length(sol.t) < 60
    @test SciMLBase.successful_retcode(sol)

    ### ROS2PR
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS2PR())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, ROS2PR())
    @test length(sol.t) < 30
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS2PR())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, ROS2PR())
    @test length(sol.t) < 30
    @test SciMLBase.successful_retcode(sol)

    ### ROS2S
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS2S())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, ROS2S())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS2S())
    @test sim.𝒪est[:final] ≈ 2 atol = testTol

    sol = solve(prob, ROS2S())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS3
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS3PR
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3PR())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3PR())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol) #length(sol.t) = 4 => Too Small??

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3PR())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3PR())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol) #length(sol.t) = 4 => Too Small??

    ### Scholz4_7
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Scholz4_7())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, Scholz4_7())
    @test length(sol.t) < 30
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Scholz4_7())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, Scholz4_7())
    @test length(sol.t) < 30
    @test SciMLBase.successful_retcode(sol)

    println("4th order Rosenbrocks")

    ### RosShamp4

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, RosShamp4())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, RosShamp4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, RosShamp4())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, RosShamp4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### Veldd4

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Veldd4())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Veldd4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Veldd4())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Veldd4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### Velds4

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Velds4())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Velds4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Velds4())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Velds4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### GRK4T

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, GRK4T())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, GRK4T())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, GRK4T())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, GRK4T())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### GRK4A
    dts = (1 / 2) .^ (7:-1:4)

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, GRK4A())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, GRK4A())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, GRK4A())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, GRK4A())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### Ros4LStab

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Ros4LStab())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Ros4LStab())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Ros4LStab())
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    sol = solve(prob, Ros4LStab())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### Rosenbrock-W Algorithms

    println("Rosenbrock-W")

    ### ROS34PW1a
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW1a(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PW1a())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW1a(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PW1a())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS34PW1b
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW1b(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PW1b())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW1b(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PW1b())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS34PW2
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW2(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PW2())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW2(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PW2())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS34PRw
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PRw())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PRw())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PRw())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS34PRw())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS3PRL
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3PRL())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3PRL())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3PRL())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3PRL())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### ROS3PRL2
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3PRL2())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3PRL2())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3PRL2())
    @test sim.𝒪est[:final] ≈ 3 atol = testTol

    sol = solve(prob, ROS3PRL2())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### RosenbrockW6S4OS
    sim = test_convergence(dts, prob, RosenbrockW6S4OS()) #test inplace
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, RosenbrockW6S4OS()) #test non-inplace
    @test sim.𝒪est[:final] ≈ 4 atol = testTol

    println("Rodas3P")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas3P(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, Rodas3P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas3P(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 3 atol = testTol
    @test sim.𝒪est[:L2] ≈ 3 atol = testTol

    sol = solve(prob, Rodas3P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sim = test_convergence(
            dts,
            prob,
            Rodas3P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test sim.𝒪est[:final] ≈ 3 atol = testTol

        sol = solve(
            prob,
            Rodas3P(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                ),
                linsolve = LinearSolve.KrylovJL()
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    ### Rodas4 Algorithms

    println("RODAS")

    dts = (1 / 2) .^ (7:-1:4)

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas4(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff()), dense_errors = true
    )
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff()))
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    if isempty(VERSION.prerelease)
        sol = solve(
            prob,
            Rodas4(
                autodiff = AutoEnzyme(
                    mode = set_runtime_activity(Enzyme.Forward), function_annotation = Enzyme.Const
                )
            )
        )
        @test length(sol.t) < 20
        @test SciMLBase.successful_retcode(sol)
    end

    sim = test_convergence(dts, prob, Rodas42(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 5.1 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas42())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(dts, prob, Rodas4P(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(dts, prob, Rodas4P2(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4P2())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas4(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    println("Rodas4 with finite diff")

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff()), dense_errors = true
    )
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff()))
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:forward))),
        dense_errors = true
    )
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:forward))))
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:complex))),
        dense_errors = true
    )
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:forward))))
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(dts, prob, Rodas42(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 5 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas42())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(dts, prob, Rodas4P(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4P())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    sim = test_convergence(dts, prob, Rodas4P2(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4P2())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    println("Rodas4P2 with finite diff")

    sim = test_convergence(
        dts, prob, Rodas4P2(autodiff = AutoFiniteDiff()), dense_errors = true
    )
    @test sim.𝒪est[:final] ≈ 4 atol = testTol
    @test sim.𝒪est[:L2] ≈ 4 atol = testTol

    sol = solve(prob, Rodas4P2(autodiff = AutoFiniteDiff()))
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    ### Rodas5
    println("Rodas5")

    prob = prob_ode_linear

    dts = (1 / 2) .^ (5:-1:2)
    sim = test_convergence(dts, prob, Rodas5(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 5 atol = testTol
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5(), dense_errors = true)
    @test sim.𝒪est[:final] ≈ 5 atol = testTol
    @test sim.𝒪est[:L2] ≈ 5 atol = testTol

    sol = solve(prob, Rodas5())
    @test length(sol.t) < 20
    @test SciMLBase.successful_retcode(sol)

end

@testset "ADTypes" begin
    for T in [
            Rosenbrock32,
            RosShamp4,
            Veldd4,
            Velds4,
            GRK4T,
            GRK4A,
            Ros4LStab,
            ROS3P,
            Rodas3,
            Rodas3P,
            Rodas4,
            Rodas42,
            Rodas4P,
            Rodas4P2,
            Rodas5,
            RosenbrockW6S4OS,
            ROS34PW1a,
            ROS34PW1b,
            ROS34PW2,
            ROS34PRw,
            ROS3PRL,
            ROS3PRL2,
            ROS2,
            ROS2PR,
            ROS2S,
            ROS3,
            ROS3PR,
            Scholz4_7,
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

        # Test with AutoFiniteDiff
        ad = AutoFiniteDiff(; fdtype = Val{:central}())
        alg = T(; autodiff = ad)
        @test alg isa RosenbrockAlgorithm
        @test alg.autodiff isa AutoFiniteDiff{Val{:central}}

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
