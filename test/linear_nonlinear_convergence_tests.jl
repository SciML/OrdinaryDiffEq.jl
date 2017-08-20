using OrdinaryDiffEq, Base.Test, DiffEqDevTools, SpecialMatrices
const μ = 1.01
f2 = (t,u) -> μ * u
f1 = DiffEqArrayOperator(μ)
(p::typeof(f1))(::Type{Val{:analytic}},t,u0) = u0.*exp.(2μ*t)

prob = SplitODEProblem(f1,f2,1/2,(0.0,1.0))
srand(100)
dts = 1./2.^(7:-1:4) #14->7 good plot
println("IIF scalar")
sim  = test_convergence(dts,prob,IIF1())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,IIF2())
@test abs(sim.𝒪est[:l2]-2) < 0.2
sim  = test_convergence(dts,prob,LawsonEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2
sim  = test_convergence(dts,prob,NorsettEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2

u0 = rand(2)
A = Strang(2)
f1 = DiffEqArrayOperator(A)
f2 = (t,u,du) -> du .= μ .* u
function (p::typeof(f1))(::Type{Val{:analytic}},t,u0)
 tmp = (A+μ*I)*t
 expm(tmp)*u0
end
prob = SplitODEProblem(f1,f2,u0,(0.0,1.0))

integrator = init(prob,NorsettEuler(),dt=1/10)
step!(integrator)
integrator.cache

dts = 1./2.^(8:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,IIF1())
@test abs(sim.𝒪est[:l2]-1) < 0.2

sim  = test_convergence(dts,prob,IIF2())
@test abs(sim.𝒪est[:l2]-2) < 0.1

sim  = test_convergence(dts,prob,LawsonEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.1

sim  = test_convergence(dts,prob,NorsettEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.1
