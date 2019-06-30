using OrdinaryDiffEq, ParameterizedFunctions, Plots, DiffEqDevTools
gr() #gr(fmt=:png)
DIAGRAMS = Dict()
f = @ode_def Orego begin
  dy1 = p1*(y2+y1*(1-p2*y1-y2))
  dy2 = (y3-(1+y1)*y2)/p1
  dy3 = p3*(y1-y3)
end p1 p2 p3

p = [77.27,8.375e-6,0.161]
prob = ODEProblem(f,[1.0,2.0,3.0],(0.0,30.0),p)
sol = solve(prob,Rodas5(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

# High Tolerance
abstols = 1 ./ 10 .^ (5:8)
reltols = 1 ./ 10 .^ (1:4);
setups = [Dict(:alg=>Rosenbrock23()),
          Dict(:alg=>Rodas3()),
          Dict(:alg=>TRBDF2())]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      save_everystep=false,appxsol=test_sol,maxiters=Int(1e5))
DIAGRAMS["OregoHighTol1"] = plot(wp)


wp = WorkPrecisionSet(prob,abstols,reltols,setups;dense = false,verbose=false,
                      appxsol=test_sol,maxiters=Int(1e5),error_estimate=:l2)
DIAGRAMS["OregoHighTol2"] = plot(wp)

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      appxsol=test_sol,maxiters=Int(1e5),error_estimate=:L2)
DIAGRAMS["OregoHighTol3"] = plot(wp)

# Low Tolerance

abstols = 1 ./ 10 .^ (7:13)
reltols = 1 ./ 10 .^ (4:10)

setups = [Dict(:alg=>GRK4A()),
          Dict(:alg=>Rodas4P()),
          Dict(:alg=>Rodas4())
]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      save_everystep=false,appxsol=test_sol,maxiters=Int(1e5))
DIAGRAMS["OregoLowTol1"] = plot(wp)

wp = WorkPrecisionSet(prob,abstols,reltols,setups;verbose=false,
                      dense=false,appxsol=test_sol,maxiters=Int(1e5),error_estimate=:l2)
DIAGRAMS["OregoLowTol2"] = plot(wp)

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      appxsol=test_sol,maxiters=Int(1e5),error_estimate=:L2)
DIAGRAMS["OregoLowTol3"] = plot(wp)

