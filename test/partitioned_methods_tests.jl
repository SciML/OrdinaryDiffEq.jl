using OrdinaryDiffEq, Base.Test, RecursiveArrayTools, DiffEqDevTools

u0 = zeros(2)
v0 = ones(2)
f1 = function (t,u,v,du)
  du .= v
end
f2 = function (t,u,v,dv)
  dv .= -u
end
function (::typeof(f2))(::Type{Val{:analytic}}, x, y0)
  u0, v0 = y0
  ArrayPartition(u0*cos(x) + v0*sin(x), -u0*sin(x) + v0*cos(x))
end

prob = DynamicalODEProblem(f1,f2,u0,v0,(0.0,5.0))

sol = solve(prob,SymplecticEuler(),dt=1/2)
sol_verlet = solve(prob,VelocityVerlet(),dt=1/100)
sol_ruth3 = solve(prob,Ruth3(),dt=1/100)

interp_time = 0:0.001:5
interp = sol(0.5)
interps = sol(interp_time)


prob = SecondOrderODEProblem(f2,u0,v0,(0.0,5.0))
(::typeof(prob.f))(::Type{Val{:analytic}},t,u0) = f2(Val{:analytic},t,u0)

sol2 = solve(prob,SymplecticEuler(),dt=1/2)
sol2_verlet = solve(prob,VelocityVerlet(),dt=1/100)
sol2_ruth3 = solve(prob,Ruth3(),dt=1/100)

sol2_verlet(0.1)

@test sol[end][1] == sol2[end][1]
@test sol_verlet[end][1] == sol2_verlet[end][1]
@test sol_ruth3[end][1] == sol2_ruth3[end][1]
@test sol[end][3] == sol2[end][3]
@test sol_verlet[end][3] == sol2_verlet[end][3]
@test sol_ruth3[end][3] == sol2_ruth3[end][3]

dts = 1.//2.^(6:-1:3)
# Symplectic Euler
sim = test_convergence(dts,prob,SymplecticEuler(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 1 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 1 rtol = 1e-1
# Verlet
sim = test_convergence(dts,prob,VelocityVerlet(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1e-1
# Test that position converges faster for Verlet
position_error = :final => [mean(sim[i].u[2].x[1] - sim[i].u_analytic[2].x[1]) for i in 1:length(sim)]
@test first(DiffEqDevTools.calc𝒪estimates(position_error).second) ≈ 3.0 rtol=1e-1

# 2nd Order Tableaus
sim = test_convergence(dts,prob,VerletLeapfrog(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1e-1
sim = test_convergence(dts,prob,PseudoVerletLeapfrog(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1e-1
sim = test_convergence(dts,prob,McAte2(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 2 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 2 rtol = 1e-1

# Ruth
sim = test_convergence(dts,prob,Ruth3(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1e-1
sim = test_convergence(dts,prob,McAte3(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1e-1
sim = test_convergence(dts,prob,CandyRoz4(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
sim = test_convergence(dts,prob,McAte4(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
sim = test_convergence(dts,prob,McAte42(),dense_errors=true)
@test_broken sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test_broken sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
sim = test_convergence(dts,prob,CalvoSanz4(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1

dts = 1.//2.^(4:-1:0)
sim = test_convergence(dts,prob,McAte5(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1

sim = test_convergence(dts,prob,Yoshida6(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4.69 rtol = 1e-1
sim = test_convergence(dts,prob,KahanLi6(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1

sim = test_convergence(dts,prob,McAte8(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1

sim = test_convergence(dts,prob,KahanLi8(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1

# Nyström method
dts = 1.//2.^(9:-1:6)
sim = test_convergence(dts,prob,Nystrom4(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
sim = test_convergence(dts,prob,Nystrom4VelocityIndependent(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
sim = test_convergence(dts,prob,IRKN3(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1e-1
sim = test_convergence(dts,prob,IRKN4(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
dts = 1.0./2.0.^(5:-1:0)
sim = test_convergence(dts,prob,Nystrom5VelocityIndependent(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 5 rtol = 1e-1

dts = 1.0./2.0.^(2:-1:-2)
sim = test_convergence(dts,prob,SofSpa10(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 10 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1e-1

f = function (t,u,du)
  du.x[1] .= u.x[2]
  du.x[2] .= -2u.x[1]
end

u = ArrayPartition((u0,v0))

prob = ODEProblem(f,u,(0.0,5.0))

sol = solve(prob,Euler(),dt=1/100)

################# Out of place symplectic

using OrdinaryDiffEq, Base.Test, RecursiveArrayTools, DiffEqDevTools

u0 = 0.0
v0 = 1.0
f1 = function (t,u,v)
  v
end
f2 = function (t,u,v)
  dv = -u
end
function (::typeof(f2))(::Type{Val{:analytic}}, x, y0)
  u0, v0 = y0
  ArrayPartition(u0*cos(x) + v0*sin(x), -u0*sin(x) + v0*cos(x))
end

prob = ODEProblem((f1,f2),(u0,v0),(0.0,5.0)) # iip wrong
@test_broken sol = solve(prob,SymplecticEuler(),dt=1/2)

# Methods need BigFloat to test convergence rate
dts = big"1.0"./big"2.0".^(6:-1:2)
prob.u0 = [big"0.0", big"0.0"], [big"1.0", big"1.0"]
sim = test_convergence(dts,prob,DPRKN6(),dense_errors=true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1e-1
@test sim.𝒪est[:L2] ≈ 6 rtol = 2e-1

