using OrdinaryDiffEq, Test
import ODEProblemLibrary: prob_ode_2Dlinear
prob = prob_ode_2Dlinear

integrator = init(prob,Euler();dt=1//2^(4))
integrator = init(prob,Tsit5();dt=1//2^(4))

solve!(integrator)

sol =solve(prob,Euler();dt=1//2^(4),maxiters=Inf)
sol =solve(prob,Tsit5();dt=1//2^(20),progress=true,adaptive=false,maxiters=Inf)

# plot(sol,plot_analytic=true)

sol =solve(prob,ExplicitRK();dt=1//2^(4))

# out-of-place 2D prob to test vec & reshape
rhs_oop(u,p,t) = u
ft0=ones(2, 2)
tspan=(0.0,30.0)
prob=ODEProblem(rhs_oop,ft0,tspan)
@test_nowarn solve(prob,Rosenbrock23())
@test_nowarn solve(prob,Rodas4())
@test_nowarn solve(prob,Kvaerno5())
