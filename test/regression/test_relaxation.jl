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

using Plots

sol1 = solve(probnum, Tsit5())
sol2 = solve(probnum, Tsit5_for_relaxation())

sol3 = solve(prob, Tsit5())
sol4 = solve(prob, Tsit5_for_relaxation())


plot(sol1, label = "Old")
plot!(sol2, label = "New")



#########################################################################
##                            Trying relaxation step

#using Optimization
#using OptimizationOptimJL
using Roots
using LinearAlgebra
using UnPack

struct Relaxation{OPT, INV}
    opt::OPT
    invariant::INV
end


function (r::Relaxation)(integrator)

    @unpack t, dt, uprev, u_propose = integrator

    # We fix here the bounds of interval where we are going to look for the relaxation
    (gamma_min, gamma_max) = apriori_bounds_dt(integrator) ./ dt
    
    S_u = u_propose - uprev

    # Minimization
    # first method tried (seems to not work)
    #=
    prob_optim = OptimizationProblem(
        (gamma,p) -> norm(r.invariant(gamma[1]*S_u .+ uprev) .- r.invariant(uprev)), 
        [dt];
        lb = [gamma_min], 
        ub = [gamma_max])
    gamma_opt = solve(prob_optim, r.opt).u[1]
    =#
    # second method
    gamma_min = 0.5
    gamma_max = 1.5
    terminate_integrator = false
    if (r.invariant(gamma_min*S_u .+ uprev) .- r.invariant(uprev)) * (r.invariant(gamma_max*S_u .+ uprev) .- r.invariant(uprev)) > 0
        gamma_opt = one(td)
        terminate_integrator = true
    else
        gamma_opt = find_zero(  gamma -> r.invariant(gamma*S_u .+ uprev) .- r.invariant(uprev),
                            (max(gamma_min,0.5), min(gamma_max,1.5)),
                            r.opt())
    end

     # Updates
    change_dt!(integrator, gamma_opt*dt)
    change_u!(integrator, uprev + gamma_opt*S_u)

    #if terminate_integrator
    #    terminate!(integrator)
    #end
end

#r = Relaxation(SAMIN(), x->x.^2)

## Tests relaxation on problem

#=

# Harmonic Oscillator
f_oscillator = (u, p, t) -> [-u[2],u[1]]
prob_oscillator = ODEProblem(
    ODEFunction(f_oscillator; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))

r_oscillator = Relaxation(AlefeldPotraShi, x-> norm(x))

sol_oscillator = solve(prob_oscillator, Tsit5_for_relaxation(); modif = r, maxiters = 5)
sol_exact = [prob_oscillator.f.analytic(prob_oscillator.u0, prob_oscillator.p, t) for t in sol_oscillator.t]
niter = length(sol_oscillator.t)

plot(sol_oscillator)
plot!(sol_oscillator.t, [sol_exact[i][1] for i in 1:niter], label = "exact u[1]", lw = 4)
plot!(sol_oscillator.t, [sol_exact[i][2] for i in 1:niter], label = "exact u[2]", lw = 4)

=#

# Non Linear Oscillator
#=
f_nonlinear_oscillator = (u, p, t) -> [-u[2]/(u[1]^2 + u[2]^2),u[1]/(u[1]^2 + u[2]^2)]
prob_nonlinear_oscillator = ODEProblem(
    ODEFunction(f_nonlinear_oscillator; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))

r_nonlinear_oscillator = Relaxation(AlefeldPotraShi, x-> norm(x))

sol_nonlinear_oscillator = solve(prob_nonlinear_oscillator, Tsit5_for_relaxation())
sol_exact = [prob_oscillator.f.analytic(prob_nonlinear_oscillator.u0, prob_nonlinear_oscillator.p, t) for t in sol_nonlinear_oscillator.t]
niter = length(sol_nonlinear_oscillator.t)


# Non Linear Pendulum
f_nonlinear_pendulum = (u, p, t) -> [-sin(u[2]), u[1]]
prob_nonlinear_pendulum = ODEProblem(
    f_nonlinear_pendulum,
    [1.0, 0.0],
    (0.0, 1.0))

r_nonlinear_pendulum = Relaxation(AlefeldPotraShi, x-> x[1]^2/2 -  cos(x[2]))

sol_nonlinear_pendulum = solve(prob_nonlinear_pendulum, Tsit5_for_relaxation())

sol_ref = solve(prob_nonlinear_pendulum, Vern9())

niter = length(sol_nonlinear_pendulum.t)
=#