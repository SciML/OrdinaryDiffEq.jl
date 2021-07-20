using OrdinaryDiffEq,DiffEqProblemLibrary, DiffEqDevTools, Test

using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear

prob = prob_ode_2Dlinear
sol =solve(prob,Rosenbrock32(),dt=1/2^4)

sol =solve(prob,ExplicitRK(tableau = constructBogakiShampine3()),dt=1/2^4)
val1 = maximum(abs.(sol.u[end] - sol.u_analytic[end]))

sol2 =solve(prob,ExplicitRK(tableau = constructDormandPrince()),dt=1/2^4)
val2 = maximum(abs.(sol2.u[end] - sol2.u_analytic[end]))

sol3 =solve(prob,ExplicitRK(tableau = constructRKF8(Float64)),dt=1/2^4)
val3 = maximum(abs.(sol3.u[end] - sol3.u_analytic[end]))

@test length(sol.t)>length(sol2.t)>=length(sol3.t)
@test max(val1,val2,val3)<2e-3

function lorenz(u,p,t)
  [10.0(u[2]-u[1])
  u[1]*(28.0-u[3]) - u[2]
  u[1]*u[2] - (8/3)*u[3]]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem{false}(lorenz,u0,tspan)
sol = solve(prob,QNDF())
@test length(sol.t) < 5000
sol = solve(prob,FBDF())
@test length(sol.t) < 6600

function lorenz(du,u,p,t)
  du[1] = 10.0(u[2]-u[1])
  du[2] = u[1]*(28.0-u[3]) - u[2]
  du[3] = u[1]*u[2] - (8/3)*u[3]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem{true}(lorenz,u0,tspan)
sol = solve(prob,QNDF())
@test length(sol.t) < 5000
sol = solve(prob,FBDF())
@test length(sol.t) < 6600
