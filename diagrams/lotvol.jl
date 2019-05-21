using OrdinaryDiffEq, DiffEqDevTools, ParameterizedFunctions

f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(f,[1.0;1.0],(0.0,10.0),p)

abstols = 1 ./ 10 .^ (6:13)
reltols = 1 ./ 10 .^ (3:10);
sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

setups = [Dict(:alg=>DP5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>Vern6())
]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;appxsol=test_sol,save_everystep=false,numruns=1000,maxiters=10000)
using Plots; gr()
DIAGRAMS["LotkaVolterra"] = plot(wp);