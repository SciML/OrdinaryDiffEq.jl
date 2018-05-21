using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase,
      DiffEqProblemLibrary, Base.Test
gc()

probArr = Vector{ODEProblem}(2)

probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear

function fixed_step_ϕstar(k)
  ∇ = Vector{typeof(k[end][1])}(3)
  ∇[1] = k[end][1]
  ∇[2] = ∇[1] - k[end-1][1]
  ∇[3] = ∇[2] - k[end-1][1] + k[end-2][1]
  return ∇
end


for i = 1:2
  prob = probArr[i]
  integrator = init(prob,VCAB3(),dt=1//256,adaptive=false)
  for i = 1:3
    step!(integrator)
  end
  @test integrator.cache.g == [1,1/2,5/12,3/8]
  step!(integrator)
  @test integrator.cache.g == [1,1/2,5/12,3/8]
  step!(integrator)
  @test integrator.cache.g == [1,1/2,5/12,3/8]
end

for i = 1:2
  prob = probArr[i]
  integrator = init(prob,VCAB3(),dt=1//256,adaptive=false)
  for i = 1:3
    step!(integrator)
  end
  @test integrator.cache.ϕstar_n == fixed_step_ϕstar(integrator.sol.k)
  step!(integrator)
  @test integrator.cache.ϕstar_n == fixed_step_ϕstar(integrator.sol.k)
  step!(integrator)
  @test integrator.cache.ϕstar_n == fixed_step_ϕstar(integrator.sol.k)
end
