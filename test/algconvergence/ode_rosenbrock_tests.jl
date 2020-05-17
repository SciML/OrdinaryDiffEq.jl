@testset "Rosenbrock Tests" begin
## Breakout these since no other test of their adaptivity

using OrdinaryDiffEq, DiffEqDevTools, Test, LinearAlgebra
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear,
                              prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear

dts = (1/2) .^ (6:-1:3)
testTol = 0.2

function _test_adaptive(alg,testTol=0.2)
    dts = (1/2) .^ (6:-1:3)
    prob = prob_ode_linear

    sim = test_convergence(dts,prob,alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    
    sol = solve(prob,alg)
    @test length(sol) < 20
    
    prob = prob_ode_2Dlinear
    
    sim = test_convergence(dts,prob,alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEq.alg_order(alg) atol=testTol
    
    sol = solve(prob,alg)
    @test length(sol) < 20
end

### Rosenbrock23()

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rosenbrock23())
@test sim.𝒪est[:final] ≈ 2 atol=testTol

sol = solve(prob,Rosenbrock23())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rosenbrock23())
@test sim.𝒪est[:final] ≈ 2 atol=testTol

sol = solve(prob,Rosenbrock23())
@test length(sol) < 20

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,Rosenbrock23(linsolve=LinSolveFactorize(qr!)))
@test sim.𝒪est[:final] ≈ 2 atol=testTol

sol = solve(prob,Rosenbrock23(linsolve=LinSolveFactorize(qr!)))
@test length(sol) < 20

_test_adaptive(Rosenbrock32())
_test_adaptive(ROS3P())
_test_adaptive(Rodas3())

println("4th order Rosenbrocks")

_test_adaptive(RosShamp4())
_test_adaptive(Veldd4())
_test_adaptive(Velds4())
_test_adaptive(GRK4T())
_test_adaptive(GRK4A(),0.3)
_test_adaptive(Ros4LStab())

### Rosenbrock-W Algorithms

println("Rosenbrock-W")
_test_adaptive(ROSWASSP3P3S1C())
_test_adaptive(ROS34PW1a())
_test_adaptive(ROS34PW1b())
_test_adaptive(ROS34PW2())
_test_adaptive(ROS34PW3())

### RosenbrockW6S4OS
sim = test_convergence(dts,prob,RosenbrockW6S4OS())#test inplace
@test sim.𝒪est[:final] ≈ 4 atol=testTol

prob = prob_ode_linear

sim = test_convergence(dts,prob,RosenbrockW6S4OS())#test non-inplace
@test sim.𝒪est[:final] ≈ 4 atol=testTol

### Rodas4 Algorithms

println("RODAS")

dts = (1/2) .^ (7:-1:4)

prob = prob_ode_linear

sim = test_convergence(dts,prob,Rodas4(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4(autodiff=false),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4(autodiff=false))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas42(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 5 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas42())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4P(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4P())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rodas4(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4())
@test length(sol) < 20

println("Rodas4 with finite diff")

sim = test_convergence(dts,prob,Rodas4(autodiff=false),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4(autodiff=false))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4(autodiff=false,
                       diff_type=Val{:forward}),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4(autodiff=false,diff_type=Val{:forward}))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4(autodiff=false,
                       diff_type=Val{:complex}),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4(autodiff=false,diff_type=Val{:complex}))
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas42(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 5 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas42())
@test length(sol) < 20

sim = test_convergence(dts,prob,Rodas4P(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 4 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas4P())
@test length(sol) < 20

### Rodas5
println("Rodas5")

prob = prob_ode_linear

dts = (1/2) .^ (7:-1:3)
sim = test_convergence(dts,prob,Rodas5(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 5 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas5())
@test length(sol) < 20

prob = prob_ode_2Dlinear

sim = test_convergence(dts,prob,Rodas5(),dense_errors=true)
@test sim.𝒪est[:final] ≈ 5 atol=testTol
@test sim.𝒪est[:L2] ≈ 4 atol=testTol

sol = solve(prob,Rodas5())
@test length(sol) < 20

prob = ODEProblem((u,p,t)->0.9u, 0.1, (0., 1.0))
@test_nowarn solve(prob, Rosenbrock23(autodiff=false))
end
