using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test
prob = deepcopy(prob_ode_2Dlinear)
prob.tspan = (1.0,0.0)

sol = solve(prob,DP5(),dt=-1/4)

if !is_windows()
  sol2 = solve(prob,dopri5(),dt=-1/4)
  @test sol.t ≈ sol2.t
end

sol = solve(prob,DP5(),dt=-1/4,tstops=[0.5])
@test sol.t == [1.0,.75,.5,0]

sol = solve(prob,RK4(),dt=-1/4)
@test sol.t == [1.00,.75,.5,.25,0]

sol = solve(prob,RK4(),dt=-1/4,tstops=[0.5,0.33])
@test sol.t ≈ [1.00,.75,.5,.33,0.08,0]
