using OrdinaryDiffEq, Base.Test, DiffEqDevTools, SpecialMatrices, DiffEqOperators
u0 = rand(2)
A = DiffEqArrayOperator(Strang(2))
function (p::typeof(A))(::Type{Val{:analytic}},t,u0)
    expm(p.A*t)*u0
end

prob = ODEProblem(A,u0,(0.0,1.0))

x = rand(2)
@test A(0.0,x) == A*x

sol = solve(prob,LinearImplicitEuler())

dts = 1./2.^(8:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,LinearImplicitEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2

# using Plots; pyplot(); plot(sim)

B = ones(2)
L = AffineDiffEqOperator{Float64}((A,),(B,),rand(2))
prob = ODEProblem(L,u0,(0.0,1.0))
sol = solve(prob,LinearImplicitEuler())

B = DiffEqArrayOperator(ones(2,2))
L = AffineDiffEqOperator{Float64}((A,B),(),rand(2))
function (p::typeof(L))(::Type{Val{:analytic}},t,u0)
    expm((p.As[1].A+p.As[2].A)*t)*u0
end

prob = ODEProblem(L,u0,(0.0,1.0))
sol = solve(prob,StrangSplitting(),dt=1/10)
using Plots; pyplot; plot(sol)
