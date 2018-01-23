using OrdinaryDiffEq, Base.Test, DiffEqDevTools, SpecialMatrices, DiffEqOperators
const μ = 1.01
f2 = (u,p,t) -> μ * u
f1 = DiffEqArrayOperator(μ)
f = SplitFunction{false}(f1,f2,nothing)
prob = SplitODEProblem(f1,f2,1/2,(0.0,1.0),func_cache=1/2)
(p::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t) = u0.*exp.(2μ*t)

srand(100)
dts = 1./2.^(7:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,GenericIIF1())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,GenericIIF2())
@test abs(sim.𝒪est[:l2]-2) < 0.2
sim  = test_convergence(dts,prob,LawsonEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,NorsettEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
@test abs(sim.𝒪est[:l2]-4) < 0.2

u0 = rand(2)
A = Strang(2)
f1 = DiffEqArrayOperator(full(A))
f2 = (du,u,p,t) -> du .= μ .* u
prob = SplitODEProblem(f1,f2,u0,(0.0,1.0))
function (p::typeof(prob.f))(::Type{Val{:analytic}},u0,p,t)
 tmp = (A+μ*I)*t
 expm(tmp)*u0
end

integrator = init(prob,NorsettEuler(),dt=1/10)
step!(integrator)
integrator.cache

dts = 1./2.^(8:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,GenericIIF1())
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim  = test_convergence(dts,prob,GenericIIF2())
@test abs(sim.𝒪est[:l2]-2) < 0.1

sim  = test_convergence(dts,prob,LawsonEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,NorsettEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,ETDRK4(),dense_errors=true)
@test abs(sim.𝒪est[:l2]-4) < 0.1
@test abs(sim.𝒪est[:L2]-4) < 0.1
