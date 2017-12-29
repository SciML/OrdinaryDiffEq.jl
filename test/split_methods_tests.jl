using OrdinaryDiffEq, DiffEqDevTools, Base.Test
testTol = 0.2

# Test that the infrustructure works

f1 = (t,u) -> 2u
f2 = (t,u) -> 2u

prob = SplitODEProblem(f1,f2,1.0,(0.0,1.0))
sol = solve(prob,SplitEuler(),dt=1/10)
sol2 = solve(prob,Euler(),dt=1/10)
@test sol2[end] == sol[end]
@test sol2(0.345) == sol(0.345)


f3 = (t,u) -> 4u
prob2 = ODEProblem(f3,1.0,(0.0,1.0))
sol3 = solve(prob2,Euler(),dt=1/10)
@test sol3[end] == sol[end]
@test sol3(0.345) == sol(0.345)

u = rand(4,2)
f1 = (t,u,du) -> du.=2u
f2 = (t,u,du) -> du.=2u
prob = SplitODEProblem(f1,f2,u,(0.0,1.0))
sol = solve(prob,SplitEuler(),dt=1/10)
sol2 = solve(prob,Euler(),dt=1/10)

@test sol2[end] == sol[end]
@test sol2(0.345) == sol(0.345)

f3 = (t,u,du) -> du.=4u
prob2 = ODEProblem(f3,u,(0.0,1.0))
sol3 = solve(prob2,Euler(),dt=1/10)

@test sol3[end] == sol[end]
@test sol3(0.345) == sol(0.345)

# Now test only the first part

f1 = (t,u) -> 2u
f2 = (t,u) -> zero(u)

prob = SplitODEProblem(f1,f2,1.0,(0.0,1.0))
function (::typeof(prob.f))(::Type{Val{:analytic}},t,u0)
    exp(2t)*u0
end

sol = solve(prob,KenCarp3())
dts = 1.//2.^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1.//2.^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

# Now test only the second part

f1 = (t,u) -> zero(u)
f2 = (t,u) -> 2u

prob = SplitODEProblem(f1,f2,1.0,(0.0,1.0))
function (::typeof(prob.f))(::Type{Val{:analytic}},t,u0)
    exp(2t)*u0
end

sol = solve(prob,KenCarp3())
dts = 1.//2.^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1.//2.^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol

# Test together

f1 = (t,u) -> u
f2 = (t,u) -> u

prob = SplitODEProblem(f1,f2,1.0,(0.0,1.0))
function (::typeof(prob.f))(::Type{Val{:analytic}},t,u0)
    exp(2t)*u0
end

sol = solve(prob,KenCarp3())
dts = 1.//2.^(12:-1:8)
sim = test_convergence(dts,prob,KenCarp3())
@test abs(sim.𝒪est[:l∞]-3) < testTol

sol = solve(prob,KenCarp4())
dts = 1.//2.^(8:-1:4)
sim = test_convergence(dts,prob,KenCarp4())
@test abs(sim.𝒪est[:l∞]-4) < testTol
