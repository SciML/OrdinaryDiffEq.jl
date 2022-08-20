using OrdinaryDiffEq, DiffEqProblemLibrary, Test

import ODEProblemLibrary: prob_ode_2Dlinear,
                                               prob_ode_large2Dlinear,
                                               prob_ode_linear,
                                               prob_ode_2Dlinear_notinplace

u0=rand(300,20).*ones(300,20)/2
prob = prob_ode_2Dlinear_notinplace
prob2 = prob_ode_2Dlinear

sol =solve(prob,Euler(),dt=1//2^(4),save_everystep=false)
sol =solve(prob2,Euler(),dt=1//2^(4),save_everystep=false)

alloc1 = @allocated sol =solve(prob,Euler(),dt=1//2^(6),save_everystep=false)
alloc2 = @allocated sol2 =solve(prob2,Euler(),dt=1//2^(6),save_everystep=false)

alloc1 = @allocated sol =solve(prob,Euler(),dt=1//2^(6),save_everystep=false)
alloc2 = @allocated sol2 =solve(prob2,Euler(),dt=1//2^(6),save_everystep=false)

@test alloc2 <= alloc1

sol = solve(prob_ode_linear,Euler(),dt=1//2^(6),save_everystep=true)
sol2 = solve(prob_ode_linear,Euler(),sol[:],sol.t,sol.k;dt=1//2^(8),save_everystep=true)

sol = solve(prob_ode_large2Dlinear,Euler(),dt=1//2^(6),save_everystep=true)
sol2 = solve(prob_ode_large2Dlinear,Euler(),sol[:],sol.t,sol.k;dt=1//2^(8),save_everystep=true)

sol = solve(prob_ode_large2Dlinear,Euler(),dt=1//2^(6),save_everystep=true)
alloc1 = @allocated sol = solve(prob_ode_large2Dlinear,Euler(),dt=1//2^(8),save_everystep=true)
alloc2 = @allocated sol2 = solve(prob_ode_large2Dlinear,Euler(),sol[:],sol.t,sol.k;dt=1//2^(8),save_everystep=true)

@test alloc2 <= alloc1

sol =solve(prob,Tsit5())
