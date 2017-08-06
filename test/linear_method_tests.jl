using DiffEqBase, OrdinaryDiffEq, Base.Test, DiffEqDevTools, SpecialMatrices
u0 = rand(2)
A = Strang(2)

f = DiffEqArrayOperator(A)
function (p::typeof(f))(::Type{Val{:analytic}},t,u0)
 expm(A*t)*u0
end

prob = ODEProblem(f,u0,(0.0,1.0))

x = rand(2)
@test f(0.0,x) == A*x

sol = solve(prob,LinearImplicitEuler())

dts = 1./2.^(8:-1:4) #14->7 good plot
sim  = test_convergence(dts,prob,LinearImplicitEuler())
@test abs(sim.𝒪est[:l2]-1) < 0.2

# using Plots; pyplot(); plot(sim)
