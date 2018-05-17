using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase,
      DiffEqProblemLibrary, Base.Test
probArr = Vector{ODEProblem}(2)
probArr[1] = prob_ode_linear

probArr[2] = prob_ode_2Dlinear
srand(100)

for i = 1:2
  prob = probArr[i]
  integrator = init(prob,VCAB3(),dt=1//256,adaptive=false)
  for i = 1:3
    step!(integrator)
  end
  g = integrator.cache.g
  @test g == [1,1/2,5/12]
  step!(integrator)
  @test g == [1,1/2,5/12]
  step!(integrator)
  @test g == [1,1/2,5/12]
end
