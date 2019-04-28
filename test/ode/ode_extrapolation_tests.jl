# Import packages
using  OrdinaryDiffEq, DiffEqDevTools, Test, Random

# Define test problems
# Note that the time span in DiffEqProblemLibrary.ODEProblemLibrary is given by
# Float64 numbers

linear = (u,p,t) -> (p*u)
linear_analytic = (u0,p,t) -> u0*exp(p*t)
prob_ode_bigfloatlinear = ODEProblem(
                          ODEFunction(linear,analytic=linear_analytic),
                          big"0.5",(big"0.0",big"1.0"),big"1.01")

f_2dlinear = (du,u,p,t) -> (@. du = p*u)
f_2dlinear_analytic = (u0,p,t) -> @. u0*exp(p*t)
prob_ode_bigfloat2Dlinear = ODEProblem(
                    ODEFunction(f_2dlinear,analytic=f_2dlinear_analytic),
                  rand(BigFloat,(4,2)),(big"0.0",big"1.0"),big"1.01")

# Prepare tests
Random.seed!(100)
problem_array = [prob_ode_bigfloatlinear,prob_ode_bigfloat2Dlinear]
dts = 1 .//2 .^(8:-1:1)

testTol = 0.2

@testset "Testing extrapolation methods" begin

# Test RichardsonEuler
@testset "Testing RichardsonEuler" begin
  for prob in problem_array
    global dts

    #  Convergence test
    for j = 1:4
      sim = test_convergence(dts,prob,AitkenNeville(j,j,j))
      @test sim.𝒪est[:final] ≈ j atol=testTol
    end

     # Regression test
    sol = solve(prob,AitkenNeville(9,1,9),reltol=1e-3)
    @test length(sol.u) < 15
    sol = solve(prob,AitkenNeville(9,1,9),reltol=1e-6)
    @test length(sol.u) < 18
  end
end # AitkenNeville

# Define the subdividing sequences
sequence_array =[:harmonic, :romberg, :bulirsch]

# Test ExtrapolationMidpointDeuflhard
@testset "Testing ExtrapolationMidpointDeuflhard" begin
  for prob in problem_array, seq in sequence_array
    global dts

    # Convergence test
    for j = 1:6
      alg = ExtrapolationMidpointDeuflhard(min_order = j,
        init_order = j, max_order=j,
        sequence_symbol = seq)
      sim = test_convergence(dts,prob,alg)
      @test sim.𝒪est[:final] ≈ 2*(alg.n_init+1) atol=testTol
    end

    # TODO: Regression test
    #...

  end
end # ExtrapolationMidpointDeuflhard

# Test ExtrapolationMidpointHairerWanner
@testset "Testing ExtrapolationMidpointHairerWanner" begin
  for prob in problem_array,
     seq in sequence_array

    # Convergence test
    for j = 1:6
      alg = ExtrapolationMidpointHairerWanner(min_order = j,
        init_order = j,
        max_order=j, sequence_symbol = seq)
      sim = test_convergence(dts,prob,alg)
      @test sim.𝒪est[:final] ≈ 2(alg.n_init+1) atol=testTol
    end

    # TODO:  Regression test
    #...

  end
end # ExtrapolationMidpointHairerWanner

end # Extrapolation methods
