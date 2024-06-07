using OrdinaryDiffEq, DiffEqDevTools,  Test

import ODEProblemLibrary: prob_ode_linear,
                          prob_ode_2Dlinear

probnum = prob_ode_linear
prob = prob_ode_2Dlinear

dts = (1 / 2) .^ (7:-1:4)
testTol = 0.2


#########################################################################
# Check that the code with the new structure without modification asked 
# by the user gives the same result
# Comparaison on Tist5
### Tsit5()

println("Tsit5")
dts = (1 / 2) .^ (7:-1:3)
sim = test_convergence(dts, probnum, Tsit5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2
sim = test_convergence(dts, prob, Tsit5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2

### Tsit5() with new structure

println("Tsit5 with relaxation")
dts = (1 / 2) .^ (7:-1:3)
sim = test_convergence(dts, probnum, Tsit5_for_relaxation())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2
sim = test_convergence(dts, prob, Tsit5_for_relaxation())  
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2              

### Compareason of interpolation

sol1 = solve(probnum, Tsit5())
sol2 = solve(probnum, Tsit5_for_relaxation())

sol3 = solve(prob, Tsit5())
sol4 = solve(prob, Tsit5_for_relaxation())

using Plots
plot(sol1, label = "Old")
plot!(sol2, label = "New")
plot(sol3, label = "Old")
plot!(sol4, label = "New")


#########################################################################
##                            Trying relaxation step

using Optimization
using OptimizationOptimJL
using LinearAlgebra
using UnPack

struct Relaxation{OPT, INV}
    opt::OPT
    invariant::INV
end


function (r::Relaxation)(integrator)

    @unpack dt, uprev, u_propose = integrator
    @unpack dtmin, dtmax, tstops = integrator.opts

    # We fix here the bounds of interval where we are going to look for the relaxation
    # and taking accound the bounds [dtmin, dtmax] and the presence of tstops
    
    # Fix of dt interval
    # should be good to have a function that gives good bound for dt taking accound dtmin explicit
    # and not 
    gamma_min = dtmin  
    gamma_max = min(dtmax, first(tstops))

    S_u = dt*(u_propose-uprev) 

    ## Minimization
    prob_optim = OptimizationProblem(
        (gamma,p) -> norm(r.invariant(gamma[1]*S_u .+ uprev) .- r.invariant(uprev)), 
        [dt];
        lb = [gamma_min], 
        ub = [gamma_max])
    gamma_opt = solve(prob_optim, r.opt).u[1]

     # Updates
    change_dt!(integrator, dt * gamma_opt)
    changed_u!(integrator, uprev + gamma_opt*S_u)
end


function (r::Relaxation)(dtmin, dtmax, dt, tstops, u_propose, uprev)

    @show gamma_min = dtmin / dt  
    @show gamma_max = min(dtmax, first(tstops)) / dt
    @show S_u = dt*(u_propose-uprev) 
    target_fun(gamma_,p) = norm(r.invariant(gamma_[1]*S_u .+ uprev) .- r.invariant(uprev))
    gamma_0 = [1.0]
    @show prob_optim = OptimizationProblem(target_fun, gamma_0; lb = [gamma_min], ub = [gamma_max])
    @show gamma__opt = solve(prob_optim, r.opt).u[1]
    # new dt
    @show dt_changed = dt * gamma__opt
    @show dt_has_changed = true
    # update u
    @show uprev + dt_changed*S_u
end

r = Relaxation(SAMIN(), x->x.^2)

#x = r(0.1, 0.2, 0.15, [1], 0.8, 0.9)


## Tests relaxation

dts = BigFloat.(1 .// 2 .^ (6:-1:2))
testTol = 0.35