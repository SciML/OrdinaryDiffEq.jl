using OrdinaryDiffEq, DiffEqBase, Base.Test
srand(100)

linear = (t,u) -> (1.01*u)
(p::typeof(linear))(::Type{Val{:analytic}},t,u0) = u0*exp(1.01*t)
prob = ODEProblem(linear,1/2,(0.0,1.0))

sol =solve(prob,Tsit5(),dt=1//2^(6),tstops=[1/2])

@test 1//2 ∈ sol.t

sol =solve(prob,RK4(),dt=1//3,tstops=[1/2])

@test sol.t == [0,1/3,1/2,1/3+1/2,1]

sol =solve(prob,RK4(),dt=1//3,tstops=[1/2],d_discontinuities=[-1/2,1/2,3/2])

@test sol.t == [0,1/3,1/2,1/3+1/2,1]

integrator = init(prob,RK4(),tstops=[1/5,1/4,1/3,1/2,3/4])

sol =solve(prob,RK4(),tstops=[1/5,1/4,1/3,1/2,3/4])

@test sol.t == [0,1/5,1/4,1/3,1/2,3/4,1]

sol =solve(prob,RK4(),tstops=[0,1/5,1/4,1/3,1/2,3/4,1])

@test sol.t == [0,1/5,1/4,1/3,1/2,3/4,1]

sol = solve(prob,RK4(),tstops=0:1//16:1)

@test sol.t == collect(0:1//16:1)

sol = solve(prob,RK4(),tstops=linspace(0,1,100))

@test sol.t == collect(linspace(0,1,100))
