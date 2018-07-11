using OrdinaryDiffEq, Test, DiffEqDevTools, DiffEqOperators, Random
const μ = 1.01
linnonlin_f2 = (u,p,t) -> μ * u
linnonlin_f1 = DiffEqArrayOperator(μ)
linnonlin_f = SplitFunction(linnonlin_f1, linnonlin_f2;
                            analytic=(u0,p,t) -> u0.*exp.(2μ*t))
prob = SplitODEProblem(linnonlin_f,1/2,(0.0,1.0))

println("Out-of-place")
srand(100)
dts = 1 ./ 2 .^ (7:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,GenericIIF1())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,GenericIIF2())
@test abs(sim.𝒪est[:l2]-2) < 0.2
sim  = test_convergence(dts,prob,LawsonEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,NorsettEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,ETDRK2())
@test abs(sim.𝒪est[:l2]-2) < 0.2
sim  = test_convergence(dts,prob,ETDRK3())
@test abs(sim.𝒪est[:l2]-3) < 0.2
sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
@test abs(sim.𝒪est[:l2]-4) < 0.2
sim  = test_convergence(dts,prob,HochOst4())
@test abs(sim.𝒪est[:l2]-4) < 0.2
sim  = test_convergence(dts,prob,ETD2())
@test abs(sim.𝒪est[:l2]-2) < 0.2

println("Inplace")
u0 = rand(2)
A = [2.0 -1.0; -1.0 2.0]
linnonlin_f1 = DiffEqArrayOperator(A)
linnonlin_f2 = (du,u,p,t) -> du .= μ .* u
linnonlin_f = SplitFunction(linnonlin_f1,linnonlin_f2;
                            analytic=(u0,p,t) -> exp((A+μ*I)*t) * u0)
prob = SplitODEProblem(linnonlin_f,u0,(0.0,1.0))

dts = 1 ./ 2 .^ (8:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,GenericIIF1())
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim  = test_convergence(dts,prob,GenericIIF2())
@test abs(sim.𝒪est[:l2]-2) < 0.1

sim  = test_convergence(dts,prob,LawsonEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,NorsettEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ETDRK2())
@test abs(sim.𝒪est[:l2]-2) < 0.1

sim  = test_convergence(dts,prob,ETDRK3())
@test abs(sim.𝒪est[:l2]-3) < 0.1

sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
@test abs(sim.𝒪est[:l2]-4) < 0.1
@test abs(sim.𝒪est[:L2]-4) < 0.1

sim  = test_convergence(dts,prob,HochOst4())
@test abs(sim.𝒪est[:l2]-4) < 0.1

sim  = test_convergence(dts,prob,ETD2())
@test abs(sim.𝒪est[:l2]-2) < 0.1
