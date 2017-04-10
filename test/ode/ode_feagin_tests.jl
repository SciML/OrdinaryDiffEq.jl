using OrdinaryDiffEq, DiffEqBase, Base.Test, DiffEqProblemLibrary, DiffEqDevTools
srand(100)
u = [0.5,0.2]
prob = prob_ode_2Dlinear

## Convergence Testing
println("Convergence Test on Linear")
dts = 1.//2.^(4:-1:2)
testTol = 1

println("Feagin RKs")
sol =solve(prob,Feagin10(),dt=dts[1])

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
(p::typeof(f))(::Type{Val{:analytic}},t,u0) -> u0*exp(linear_bigα*t)
prob_ode_bigfloat2Dlinear = ODEProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,(0.0,1.0))

prob = prob_ode_bigfloat2Dlinear

sim = test_convergence(dts,prob,Feagin10())
@test abs(sim.𝒪est[:final]-8) < testTol #Lowered due to low test dt

sim = test_convergence(dts,prob,Feagin12())
@test abs(sim.𝒪est[:final]-12) < testTol

sim = test_convergence(dts,prob,Feagin14())
@test abs(sim.𝒪est[:final]-15) < testTol #Upped to 15 for test

f = (t,u) -> (linear_bigα*u)
(p::typeof(f))(::Type{Val{:analytic}},t,u0) -> u0*exp(linear_bigα*t)
prob_ode_bigfloatlinear = ODEProblem(f,parse(BigFloat,"0.5"),(0.0,1.0))
prob = prob_ode_bigfloatlinear

dts = 1.//2.^(6:-1:3)
sim = test_convergence(dts,prob,Feagin10())
@test abs(sim.𝒪est[:final]-10) < testTol

dts = 1.//2.^(4:-1:2)
sim = test_convergence(dts,prob,Feagin12())
@test abs(sim.𝒪est[:final]-12) < testTol

sim = test_convergence(dts,prob,Feagin14())
@test abs(sim.𝒪est[:final]-15) < testTol #Upped to 15 for test

prob = prob_ode_bigfloat2Dlinear

#compile
sol =solve(prob,Feagin10(),dt=dts[1])
sol =solve(prob,Feagin12(),dt=dts[1])
sol =solve(prob,Feagin14(),dt=dts[1])

#test
@time sol =solve(prob,Feagin10(),dt=dts[1])
@time sol =solve(prob,Feagin12(),dt=dts[1])
@time sol =solve(prob,Feagin14(),dt=dts[1])
