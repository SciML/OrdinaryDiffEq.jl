using OrdinaryDiffEq, Test

# in terms of the voltage across all three elements
rlc1!(v′,v,(R,L,C),t) = -(v′/R + v/L)/C
identity_f(v,u,p,t) = v # needed to form second order dynamical ODE

setup_rlc(R,L,C;v_init=0.0,v′_init=0.0,tspan=(0.0,50.0)) =
    DynamicalODEProblem{false}(rlc1!,identity_f,v′_init,v_init,tspan,(R,L,C))

# simulate voltage impulse
R,L,C = 10, 0.3, 2

prob = setup_rlc(R,L,C,v_init=2.0)

res1 = solve(prob,Vern8(),dt=1/10,saveat=1/10)
res3 = solve(prob,CalvoSanz4(),dt=1/10,saveat=1/10)

sol = solve(prob,CalvoSanz4(),dt=1/10)
@test sol(0.32) isa OrdinaryDiffEq.ArrayPartition
@test sol(0.32, Val{1}) isa OrdinaryDiffEq.ArrayPartition
@test sol(0.32, Val{2}) isa OrdinaryDiffEq.ArrayPartition
@test sol(0.32, Val{3}) isa OrdinaryDiffEq.ArrayPartition
