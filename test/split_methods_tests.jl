using OrdinaryDiffEq, DiffEqDevTools, Test, Random
testTol = 0.2

srand(100)

# Test that the infrustructure works

f1 = (u,p,t) -> 2u
f2 = (u,p,t) -> 2u

prob = SplitODEProblem(f1,f2,1.0,(0.0,1.0))
sol = solve(prob,SplitEuler(),dt=1/10)
sol2 = solve(prob,Euler(),dt=1/10)
@test sol2[end] == sol[end]
@test sol2(0.345) == sol(0.345)


f3 = (u,p,t) -> 4u
prob2 = ODEProblem(f3,1.0,(0.0,1.0))
sol3 = solve(prob2,Euler(),dt=1/10)
@test sol3[end] == sol[end]
@test sol3(0.345) == sol(0.345)

u = rand(4,2)
f1 = (du,u,p,t) -> du.=2u
f2 = (du,u,p,t) -> du.=2u
prob = SplitODEProblem(f1,f2,u,(0.0,1.0))
sol = solve(prob,SplitEuler(),dt=1/10)
sol2 = solve(prob,Euler(),dt=1/10)

@test sol2[end] == sol[end]
@test sol2(0.345) == sol(0.345)

f3 = (du,u,p,t) -> du.=4u
prob2 = ODEProblem(f3,u,(0.0,1.0))
sol3 = solve(prob2,Euler(),dt=1/10)

@test sol3[end] == sol[end]
@test sol3(0.345) == sol(0.345)

# Now test only the first part

f1 = (u,p,t) -> 2u
f2 = (u,p,t) -> zero(u)

ff_split = SplitFunction(f1, f2; analytic=(u0,p,t)->exp(2t)*u0)
prob = SplitODEProblem(ff_split,1.0,(0.0,1.0))

sol = solve(prob,KenCarp3())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp3(nonlinsolve=NLFunctional()))
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

sol = solve(prob,KenCarp5())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp5())
@test abs(sim.𝒪est[:l∞]-5) < testTol

# IMEXEuler
dts = 1 .//2 .^(8:-1:4)
sim1 = test_convergence(dts,prob,IMEXEuler())
@test abs(sim1.𝒪est[:l∞]-1) < testTol

# CNAB2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNAB2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# CNLF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNLF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF3
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF3())
@test_broken abs(sim.𝒪est[:l∞]-3) < testTol

# SBDF4
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF4())
@test_broken abs(sim.𝒪est[:l∞]-4) < testTol

# Now test only the second part

f1 = (u,p,t) -> zero(u)
f2 = (u,p,t) -> 2u

ff_split2 = SplitFunction(f1, f2; analytic=(u0,p,t)->exp(2t)*u0)
prob = SplitODEProblem(ff_split2,1.0,(0.0,1.0))

sol = solve(prob,KenCarp3())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

sol = solve(prob,KenCarp5())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp5())
@test abs(sim.𝒪est[:l∞]-5) < testTol

# IMEXEuler
dts = 1 .//2 .^(8:-1:4)
sim2 = test_convergence(dts,prob,IMEXEuler())
@test abs(sim2.𝒪est[:l∞]-1) < testTol

# CNAB2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNAB2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# CNLF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNLF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF3
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF3())
@test_broken abs(sim.𝒪est[:l∞]-3) < testTol

# SBDF4
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF4())
@test_broken abs(sim.𝒪est[:l∞]-4) < testTol

# Test together

f1 = (u,p,t) -> u
f2 = (u,p,t) -> 2u

ff_split3 = SplitFunction(f1, f2; analytic=(u0,p,t)->exp(3t)*u0)
prob = SplitODEProblem(ff_split3,1.0,(0.0,1.0))

sol = solve(prob,KenCarp3())
dts = 1 .//2 .^(12:-1:8)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

sol = solve(prob,KenCarp5())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp5())
@test abs(sim.𝒪est[:l∞]-5) < testTol

# IMEXEuler
dts = 1 .//2 .^(8:-1:4)
sim3 = test_convergence(dts,prob,IMEXEuler())
@test abs(sim3.𝒪est[:l∞]-1) < testTol

# CNAB2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNAB2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# CNLF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNLF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF3
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF3())
@test_broken abs(sim.𝒪est[:l∞]-3) < testTol

# SBDF4
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF4())
@test_broken abs(sim.𝒪est[:l∞]-4) < testTol

# Now test only the first part

f1 = (du,u,p,t) -> du .= 2u
f2 = (du,u,p,t) -> du .= 0.0

ff_split4 = SplitFunction(f1, f2; analytic=(u0,p,t)->exp(2t)*u0)
prob = SplitODEProblem(ff_split4,rand(4,2),(0.0,1.0))

sol = solve(prob,KenCarp3())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

sol = solve(prob,KenCarp5())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp5())
@test abs(sim.𝒪est[:l∞]-5) < testTol

# IMEXEuler
dts = 1 .//2 .^(8:-1:4)
sim1 = test_convergence(dts,prob,IMEXEuler())
@test abs(sim1.𝒪est[:l∞]-1) < testTol

# CNAB2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNAB2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# CNLF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNLF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF3
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF3())
@test_broken abs(sim.𝒪est[:l∞]-3) < testTol

# SBDF4
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF4())
@test_broken abs(sim.𝒪est[:l∞]-4) < testTol

# Now test only the second part

f1 = (du,u,p,t) -> du.= 0.0
f2 = (du,u,p,t) -> du.= 2u

ff_split5 = SplitFunction(f1, f2; analytic=(u0,p,t)->exp(2t)*u0)
prob = SplitODEProblem(ff_split5,rand(4,2),(0.0,1.0))

sol = solve(prob,KenCarp3())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

sol = solve(prob,KenCarp5())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp5())
@test abs(sim.𝒪est[:l∞]-5) < testTol

# IMEXEuler
dts = 1 .//2 .^(8:-1:4)
sim2 = test_convergence(dts,prob,IMEXEuler())
@test abs(sim2.𝒪est[:l∞]-1) < testTol

# CNAB2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNAB2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# CNLF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNLF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF3
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF3())
@test_broken abs(sim.𝒪est[:l∞]-3) < testTol

# SBDF4
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF4())
@test_broken abs(sim.𝒪est[:l∞]-4) < testTol

# Test together

f1 = (du,u,p,t) -> du .= u
f2 = (du,u,p,t) -> du .= 2u

ff_split6 = SplitFunction(f1, f2; analytic=(u0,p,t)->exp(3t)*u0)
prob = SplitODEProblem(ff_split6,rand(4,2),(0.0,1.0))

sol = solve(prob,KenCarp3())
dts = 1 .//2 .^(12:-1:8)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

sol = solve(prob,KenCarp5())
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp5())
@test abs(sim.𝒪est[:l∞]-5) < testTol

# IMEXEuler
dts = 1 .//2 .^(8:-1:4)
sim3 = test_convergence(dts,prob,IMEXEuler())
@test abs(sim3.𝒪est[:l∞]-1) < testTol

# CNAB2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNAB2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# CNLF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,CNLF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF2
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF2())
@test abs(sim.𝒪est[:l∞]-2) < testTol

# SBDF3
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF3())
@test_broken abs(sim.𝒪est[:l∞]-3) < testTol

# SBDF4
dts = 1 .//2 .^(8:-1:4)
sim = test_convergence(dts,prob,SBDF4())
@test_broken abs(sim.𝒪est[:l∞]-4) < testTol
