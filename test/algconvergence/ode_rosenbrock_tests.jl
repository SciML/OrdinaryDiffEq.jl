using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra, LinearSolve
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear,
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

    ### RosenbrockW6S4OS
    sim = test_convergence(dts, prob, RosenbrockW6S4OS())#test inplace
    @test sim.𝒪est[:final]≈4 atol=testTol

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, RosenbrockW6S4OS())#test non-inplace
    @test sim.𝒪est[:final]≈4 atol=testTol

    ### Rodas4 Algorithms

    println("RODAS")

    dts = (1 / 2) .^ (7:-1:4)

    prob = prob_ode_linear

    sim = test_convergence(dts, prob, Rodas4(), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4())
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4(autodiff = false), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = false))
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

    sim = test_convergence(dts, prob, Rodas4(autodiff = false), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = false))
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4(autodiff = false,
                                             diff_type = Val{:forward}),
                           dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = false, diff_type = Val{:forward}))
    @test length(sol) < 20

    sim = test_convergence(dts, prob, Rodas4(autodiff = false,
                                             diff_type = Val{:complex}),
                           dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4(autodiff = false, diff_type = Val{:complex}))
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

    sim = test_convergence(dts, prob, Rodas4P2(autodiff = false), dense_errors = true)
    @test sim.𝒪est[:final]≈4 atol=testTol
    @test sim.𝒪est[:L2]≈4 atol=testTol

    sol = solve(prob, Rodas4P2(autodiff = false))
    @test length(sol) < 20

    ### Rodas5
    println("Rodas5")

    prob = prob_ode_linear

    dts = (1 / 2) .^ (6:-1:3)
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

    dts = (1 / 2) .^ (6:-1:3)
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

    prob = ODEProblem((u, p, t) -> 0.9u, 0.1, (0.0, 1.0))
    @test_nowarn solve(prob, Rosenbrock23(autodiff = false))
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
