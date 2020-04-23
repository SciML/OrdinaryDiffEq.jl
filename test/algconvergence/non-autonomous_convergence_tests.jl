using OrdinaryDiffEq, DiffEqDevTools, Test

function nonauto1(u, p, t)
  x, _ = u
  [
    t*x
    0
  ]
end

function nonauto2(u, p, t)
  _, y = u
  [
   y,
   t*y
  ]
end

function analytic(u0,p,t)
  x0, y0 = u0
  et = exp(t^2/2)
  [
   et*(x0 + t*y0),
   et*y0
  ]
end

u0 = [1.1, 2.2]
tspan = (0.0,1.0)
prob1 = ODEProblem(ODEFunction{true}((du,u,p,t)->du .= nonauto1(u,p,t) .+ nonauto2(u,p,t),
                              analytic=analytic),
                              u0,tspan)
prob2 = ODEProblem(ODEFunction{false}((u,p,t)-> nonauto1(u,p,t) .+ nonauto2(u,p,t),
                              analytic=analytic),
                              u0,tspan)
prob3 = SplitODEProblem(SplitFunction{true}((du,u,p,t)->du .= nonauto1(u,p,t),
                                       (du,u,p,t)->du .= nonauto2(u,p,t),
                              analytic=analytic),
                              u0,tspan)
prob4 = SplitODEProblem(SplitFunction{false}(nonauto1,
                                        nonauto2,
                              analytic=analytic),
                              u0,tspan)

testTol = 0.2

kwargs = (reltol=1e-16, abstol=1e-16)
for prob in [prob1, prob2, prob3, prob4]
  dts = 1 .//2 .^(12:-1:8)
  sim = test_convergence(dts,prob,KenCarp3(); kwargs...)
  @test sim.𝒪est[:l∞] ≈ 3 atol=testTol
  dts = 1 .//2 .^(9:-1:6)
  sim = test_convergence(dts,prob,CFNLIRK3(); kwargs...)
  @test sim.𝒪est[:l∞] ≈ 3 atol=testTol
  sim = test_convergence(dts,prob,CFNLIRK4(); kwargs...)
  @test sim.𝒪est[:l∞] ≈ 4 atol=testTol
  sim = test_convergence(dts,prob,KenCarp4(); kwargs...)
  @test sim.𝒪est[:l∞] ≈ 4 atol=testTol
  dts = 1 .//2 .^(7:-1:4)
  sim = test_convergence(dts,prob,KenCarp5(); kwargs...)
  @test sim.𝒪est[:l∞] ≈ 5 atol=testTol
end
for prob in [prob1, prob2]
  dts = 1 .//2 .^(7:-1:4)
  sim = test_convergence(dts,prob,ESDIRK54I8L2SA(); kwargs...)
  @test sim.𝒪est[:l∞] ≈ 5 atol=testTol
end

#=
prob = prob2
dts = 1 .//2 .^(9:-1:6)
sim = test_convergence(dts,prob,KenCarp4())
sim.𝒪est[:l∞]

using Plots
plot(sim)
=#
