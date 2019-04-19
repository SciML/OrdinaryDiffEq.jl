# This definitely needs cleaning

# While developing use Revise
# TODO: Delete this later
using Revise

# Import packages
using  OrdinaryDiffEq, DiffEqDevTools, Test, Random

# Import test problems
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear,
  prob_ode_bigfloat2Dlinear, prob_ode_2Dlinear_notinplace

# Prepare tests
Random.seed!(100)
problem_array = [prob_ode_linear, prob_ode_2Dlinear,
  prob_ode_bigfloat2Dlinear,
  prob_ode_2Dlinear_notinplace]
dts = 1 .//2 .^(8:-1:4)

testTol = 0.2

# probArr[1] = prob_ode_linear
# probArr[2] = prob_ode_2Dlinear
@testset "Testing extrapolation methods" begin

# Test RichardsonEuler
@testset "Testing RichardsonEuler" begin
  for prob in problem_array
    global dts

    #  Convergence test
    for j = 1:4
      sim = test_convergence(dts,prob,RichardsonEuler(j,j,j))
      @test sim.𝒪est[:final] ≈ j atol=testTol
    end

     # Regression test
    sol = solve(prob,RichardsonEuler(9,1,9),reltol=1e-3)
    @test length(sol.u) < 15
    sol = solve(prob,RichardsonEuler(9,1,9),reltol=1e-6)
    @test length(sol.u) < 18
  end
end # RichardsonEuler

# Test of ExtrapolationMidpoint*() algorithms

# Define the subdividing sequences
sequence_array =[:harmonic, :romberg, :bulirsch]

# Test ExtrapolationMidpointDeuflhard
@testset "Testing ExtrapolationMidpointDeuflhard" begin
  for prob in problem_array, seq in sequence_array
    global dts

    # Convergence test
    for j = 1:4
      alg = ExtrapolationMidpointDeuflhard(min_extrapolation_order = j,
        init_extrapolation_order = j, max_extrapolation_order=j,
        sequence_symbol = seq)
      sim = test_convergence(dts,prob,alg)
      @test sim.𝒪est[:final] ≈ 2*(j+1) atol=testTol
    end

    # TODO: Regression test
    #...

  end
end # ExtrapolationMidpointDeuflhard

# Test ExtrapolationMidpointHairerWanner
@testset "Testing ExtrapolationMidpointHairerWanner" begin
  for prob in problem_array, seq in sequence_array
    global dts

    # Convergence test
    for j = 1:4
      alg = ExtrapolationMidpointHairerWanner(min_extrapolation_order = j,
        init_extrapolation_order = j,
        max_extrapolation_order=j, sequence_symbol = seq)
      sim = test_convergence(dts,prob,alg)
      @test sim.𝒪est[:final] ≈ 2(j+1) atol=testTol
    end

    # TODO:  Regression test
    #...

  end
end # ExtrapolationMidpointHairerWanner

end # Extrapolation methods
