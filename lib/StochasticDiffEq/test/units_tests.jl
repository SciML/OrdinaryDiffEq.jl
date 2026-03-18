# Stochastic needs ΔW = s^(1/2).
#=
β = 0.6
σ = (t,y) -> β*y/(4.0s)
u0 = 1.5Newton
prob = SDEProblem(f,σ,u0)

sol =solve(prob::SDEProblem,[0,1],dt=(1/2^4)Second,save_everystep=true,alg=:EM)
sol =solve(prob::SDEProblem,[0,1],dt=(1/2^4)Second,save_everystep=true,alg=:SRIW1)

TEST_PLOT && plot(sol)

u0 = [1.5Newton 2.0Newton
      3.0Newton 1.0Newton]

prob = SDEProblem(f,σ,u0)

sol =solve(prob::SDEProblem,[0,1],dt=(1/2^4)Second,save_everystep=true,alg=:EM)
sol =solve(prob::SDEProblem,[0,1],dt=(1/2^4)Second,save_everystep=true,alg=:SRIW1)
=#
