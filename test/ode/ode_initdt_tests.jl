using OrdinaryDiffEq,DiffEqProblemLibrary, DiffEqDevTools, Base.Test

prob = prob_ode_linear
println("Solve and Plot")
sol =solve(prob,Rosenbrock32())
dt₀ = sol.t[2]

prob = prob_ode_2Dlinear

## Solve and plot
println("Solve and Plot")
tab = constructBogakiShampine3()
sol =solve(prob,ExplicitRK(),tableau=tab)
dt₀ = sol.t[2]

@test  1e-7 < dt₀ < .1
@test_throws ErrorException sol = solve(prob,Euler())
#dt₀ = sol.t[2]

tab = constructDormandPrince8_64bit()
sol3 =solve(prob,ExplicitRK(),tableau=tab)
dt₀ = sol3.t[2]

@test 1e-7 < dt₀ < .3

prob = prob_ode_linear
sol = solve(prob,DP5())

if !is_windows()
  using ODEInterfaceDiffEq, Base.Test
  sol2 = solve(prob,dopri5())
  @test sol.t[2] ≈ sol2.t[2]
end


prob = prob_ode_2Dlinear
sol = solve(prob,DP5(),internalnorm=(u)->sqrt(sum(abs2,u)))

if !is_windows()
  # Change the norm due to error in dopri5.f
  sol2 = solve(prob,dopri5())
  @test sol.t[2] ≈ sol2.t[2]
end

prob = deepcopy(prob_ode_linear)
prob.tspan = (1.0,0.0)
sol = solve(prob,DP5())

if !is_windows()
  sol2 = solve(prob,dopri5())
  @test sol.t[2] ≈ sol2.t[2]
end

prob = deepcopy(prob_ode_2Dlinear)
prob.tspan = (1.0,0.0)
sol = solve(prob,DP5(),internalnorm=(u)->sqrt(sum(abs2,u)))

if !is_windows()
  # Change the norm due to error in dopri5.f
  sol2 = solve(prob,dopri5())
  @test sol.t[2] ≈ sol2.t[2]
end
