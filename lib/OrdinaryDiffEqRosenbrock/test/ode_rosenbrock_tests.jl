using OrdinaryDiffEqRosenbrock, DiffEqDevTools, Test, LinearAlgebra, LinearSolve, ADTypes
import ODEProblemLibrary: prob_ode_linear,
                          prob_ode_2Dlinear,
                          prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear
import LinearSolve

@testset "Rosenbrock Tests" begin
    ## Breakout these since no other test of their adaptivity

    dts = (1 / 2) .^ (6:-1:3)
    testTol = 0.2

    ### Rosenbrock23()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rosenbrock23())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, Rosenbrock23())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock23())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, Rosenbrock23())
    @test length(sol) < 20

    prob = prob_ode_bigfloat2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock23(linsolve = QRFactorization()))
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, Rosenbrock23(linsolve = QRFactorization()))
    @test length(sol) < 20

    ### Rosenbrock32()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rosenbrock32())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, Rosenbrock32())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rosenbrock32())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, Rosenbrock32())
    @test length(sol) < 20

    ### ROS3P()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3P())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3P())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3P())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3P())
    @test length(sol) < 20

    ### Rodas3()

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas3())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, Rodas3())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas3())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, Rodas3())
    @test length(sol) < 20

    ### ROS2
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS2())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, ROS2())
    @test length(sol) < 61

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS2())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, ROS2PR())
    @test length(sol) < 60

    ### ROS2PR
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS2PR())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, ROS2PR())
    @test length(sol) < 30

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS2PR())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, ROS2PR())
    @test length(sol) < 30

    ### ROS2S
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS2S())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, ROS2S())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS2S())
    @test sim.𝒪est[:final]≈2 atol=testTol

    sol = solve(prob, ROS2S())
    @test length(sol) < 20

    ### ROS3
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3())
    @test length(sol) < 20

    ### ROS3PR
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3PR())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3PR())
    @test length(sol) < 20 #length(sol) = 4 => Too Small??

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3PR())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3PR())
    @test length(sol) < 20 #length(sol) = 4 => Too Small??

    ### Scholz4_7
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Scholz4_7())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, Scholz4_7())
    @test length(sol) < 30

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Scholz4_7())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, Scholz4_7())
    @test length(sol) < 30

    println("4th order Rosenbrocks")

    ### RosShamp4

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, RosShamp4())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, RosShamp4())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, RosShamp4())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, RosShamp4())
    @test length(sol) < 20

    ### Veldd4

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Veldd4())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, Veldd4())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Veldd4())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, Veldd4())
    @test length(sol) < 20

    ### Velds4

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Velds4())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, Velds4())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Velds4())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, Velds4())
    @test length(sol) < 20

    ### GRK4T

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, GRK4T())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, GRK4T())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, GRK4T())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, GRK4T())
    @test length(sol) < 20

    ### GRK4A
    dts = (1 / 2) .^ (7:-1:4)

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, GRK4A())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, GRK4A())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, GRK4A())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, GRK4A())
    @test length(sol) < 20

    ### Ros4LStab

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Ros4LStab())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, Ros4LStab())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Ros4LStab())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, Ros4LStab())
    @test length(sol) < 20

    ### Rosenbrock-W Algorithms

    println("Rosenbrock-W")

    ### ROS34PW1a
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW1a())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PW1a())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW1a())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PW1a())
    @test length(sol) < 20

    ### ROS34PW1b
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW1b())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PW1b())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW1b())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PW1b())
    @test length(sol) < 20

    ### ROS34PW2
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW2())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PW2())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW2())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PW2())
    @test length(sol) < 20

    ### ROS34PW3
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PW3())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, ROS34PW3())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PW3())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, ROS34PW3())
    @test length(sol) < 20

    ### ROS34PRw
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS34PRw())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PRw())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS34PRw())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS34PRw())
    @test length(sol) < 20

    ### ROS3PRL
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3PRL())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3PRL())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3PRL())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3PRL())
    @test length(sol) < 20

    ### ROS3PRL2
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROS3PRL2())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3PRL2())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROS3PRL2())
    @test sim.𝒪est[:final]≈3 atol=testTol

    sol = solve(prob, ROS3PRL2())
    @test length(sol) < 20

    ### ROK4a
    prob = prob_ode_linear

    sim = test_convergence(dts, prob, ROK4a())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, ROK4a())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, ROK4a())
    @test sim.𝒪est[:final]≈4 atol=testTol

    sol = solve(prob, ROK4a())
    @test length(sol) < 20

    ### RosenbrockW6S4OS
    sim = test_convergence(dts, prob, RosenbrockW6S4OS())#test inplace
    @test sim.𝒪est[:final]≈4 atol=testTol

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, RosenbrockW6S4OS())#test non-inplace
    @test sim.𝒪est[:final]≈4 atol=testTol

    ### Rodas23W, Rodas3P

    println("Rodas23W")

    prob = prob_ode_linear

    dts = (1 / 2) .^ (6:-1:3)
    sim = test_convergence(dts, prob, Rodas23W(), dense_errors = true)
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:L2]≈2 atol=testTol

    sol = solve(prob, Rodas23W())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas23W(), dense_errors = true)
    @test sim.𝒪est[:final]≈2 atol=testTol
    @test sim.𝒪est[:L2]≈2 atol=testTol

    sol = solve(prob, Rodas23W())
    @test length(sol) < 20

    println("Rodas3P")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas3P(), dense_errors = true)
    @test sim.𝒪est[:final]≈3 atol=testTol
    @test sim.𝒪est[:L2]≈3 atol=testTol

    sol = solve(prob, Rodas3P())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas3P(), dense_errors = true)
    @test sim.𝒪est[:final]≈3 atol=testTol
    @test sim.𝒪est[:L2]≈3 atol=testTol

    sol = solve(prob, Rodas3P())
    @test length(sol) < 20

    ### Rodas4 Algorithms

    println("RODAS")

    dts = (1 / 2) .^ (7:-1:4)

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas4(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4())
    @test length(sol) < 20

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff()), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff()))
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas42(), dense_errors = true)
    @test sim.𝒪est[:final]≈5.1 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas42())
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4P(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4P())
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4P2(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4P2())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas4(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4())
    @test length(sol) < 20

    println("Rodas4 with finite diff")

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff()), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff()))
    @test length(sol) < 20

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:forward))),
        dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:forward))))
    @test length(sol) < 20

    sim = test_convergence(
        dts, prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:complex))),
        dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = AutoFiniteDiff(fdtype = Val(:forward))))
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas42(), dense_errors = true)
    @test sim.𝒪est[:final]≈5 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas42())
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4P(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4P())
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4P2(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4P2())
    @test length(sol) < 20

    println("Rodas4P2 with finite diff")

    sim = test_convergence(
        dts, prob, Rodas4P2(autodiff = AutoFiniteDiff()), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4P2(autodiff = AutoFiniteDiff()))
    @test length(sol) < 20

    ### Rodas5
    println("Rodas5")

    prob = prob_ode_linear

    dts = (1 / 2) .^ (5:-1:2)
    sim = test_convergence(dts, prob, Rodas5(), dense_errors = true)
    @test sim.𝒪est[:final]≈5 atol=testTol
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5(), dense_errors = true)
    @test sim.𝒪est[:final]≈5 atol=testTol
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5())
    @test length(sol) < 20

    println("Rodas5P")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5P())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5P(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5P())
    @test length(sol) < 20

    println("Rodas5Pe")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5Pe(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5Pe())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5Pe(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5Pe())
    @test length(sol) < 20

    println("Rodas5Pr")

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas5Pr(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5Pr())
    @test length(sol) < 20

    prob = prob_ode_2Dlinear

    sim = test_convergence(dts, prob, Rodas5Pr(), dense_errors = true)
    #@test sim.𝒪est[:final]≈5 atol=testTol #-- observed order > 6
    @test sim.𝒪est[:L2]≈5 atol=testTol

    sol = solve(prob, Rodas5Pr())
    @test length(sol) < 20

    prob = ODEProblem((u, p, t) -> 0.9u, 0.1, (0.0, 1.0))
    @test_nowarn solve(prob, Rosenbrock23(autodiff = AutoFiniteDiff()))
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
    sim = test_convergence(dts, prob, Rodas3(linsolve = LinearSolve.KrylovJL()))
    @test sim.𝒪est[:final]≈3 atol=testTol
end

@testset "ADTypes" begin
    for T in [
        Rosenbrock23,
        Rosenbrock32,
        RosShamp4,
        Veldd4,
        Velds4,
        GRK4T,
        GRK4A,
        Ros4LStab,
        ROS3P,
        Rodas3,
        Rodas23W,
        Rodas3P,
        Rodas4,
        Rodas42,
        Rodas4P,
        Rodas4P2,
        Rodas5,
        Rodas5P,
        Rodas5Pe,
        Rodas5Pr,
        RosenbrockW6S4OS,
        ROS34PW1a,
        ROS34PW1b,
        ROS34PW2,
        ROS34PW3,
        ROS34PRw,
        ROS3PRL,
        ROS3PRL2,
        ROK4a,
        ROS2,
        ROS2PR,
        ROS2S,
        ROS3,
        ROS3PR,
        Scholz4_7
    ]
        RosenbrockAlgorithm = if T <:
                                 OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAlgorithm
            OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAlgorithm
        else
            OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
        end

        ad = AutoForwardDiff(; chunksize = 3)
        alg = @test_logs @inferred(T(; autodiff = ad))
        @test alg isa RosenbrockAlgorithm{3, typeof(ad), Val{:forward}()}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.alg_autodiff(alg) === ad
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.get_chunksize(alg) === Val{3}()

        alg = @test_logs (:warn, r"The `chunk_size` keyword is deprecated") match_mode=:any @inferred(T(;
            autodiff = ad, chunk_size = Val{4}()))
        @test alg isa RosenbrockAlgorithm{4, <:AutoForwardDiff{4}, Val{:forward}()}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.alg_autodiff(alg) isa
              AutoForwardDiff{4}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.get_chunksize(alg) === Val{4}()

        ad = AutoFiniteDiff(; fdtype = Val{:central}())
        alg = @test_logs @inferred(T(; autodiff = ad))
        @test alg isa
              RosenbrockAlgorithm{0, <:AutoFiniteDiff{Val{:central}}, Val{:central}()}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.alg_autodiff(alg) === ad
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.get_chunksize(alg) === Val{0}()

        alg = @test_logs (:warn, r"The `diff_type` keyword is deprecated") match_mode=:any @inferred(T(;
            autodiff = ad, diff_type = Val{:complex}()))
        @test alg isa
              RosenbrockAlgorithm{0, <:AutoFiniteDiff{Val{:complex}}, Val{:complex}()}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.alg_autodiff(alg) isa
              AutoFiniteDiff{Val{:complex}}
        @test OrdinaryDiffEqRosenbrock.OrdinaryDiffEqCore.get_chunksize(alg) === Val{0}()

        # issue #2613
        f(u, _, _) = -u
        prob = ODEProblem(f, [1.0, 0.0], (0.0, 1.0))
        alg = T(; autodiff = AutoForwardDiff(; chunksize = 1))
        sol = if alg isa OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
            @inferred(solve(prob, alg))
        else
            @inferred(solve(prob, alg; dt = 0.1))
        end
        @test sol.alg === alg
        alg = T(; autodiff = AutoFiniteDiff(; fdtype = Val(:central)))
        sol = if alg isa OrdinaryDiffEqRosenbrock.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
            @inferred(solve(prob, alg))
        else
            @inferred(solve(prob, alg; dt = 0.1))
        end
        @test sol.alg === alg
    end
end
