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
    #(gamma_min, gamma_max) = apriori_bounds_dt(integrator) ./ dt
    
    S_u = u_propose - uprev

    # Find relaxation paramter gamma
    gamma_min = 0.5
    gamma_max = 1.5

    if (r.invariant(gamma_min*S_u .+ uprev) .- r.invariant(uprev)) * (r.invariant(gamma_max*S_u .+ uprev) .- r.invariant(uprev)) ≤ 0
        gamma = find_zero(gamma -> r.invariant(gamma*S_u .+ uprev) .- r.invariant(uprev),
                            (gamma_min, gamma_max),
                            r.opt())

        change_dt!(integrator, gamma*dt)
        change_u!(integrator, uprev + gamma*S_u)
    end

    # println("######################  Print of dt time")
    # @show integrator.dt
    # @show integrator.dt_changed
    # @show integrator.opts.dtmin
    # @show integrator.opts.dtmax
    # if has_tstop(integrator)
    #     @show first_tstop(integrator) - integrator.t
    # end

end


########################################################
# TEST  1 : Harmonic Oscillator
printstyled("Harmonic Oscillator\n"; bold = true)

f_oscillator = (u, p, t) -> [-u[2],u[1]]
prob_oscillator = ODEProblem(
    ODEFunction(f_oscillator; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))
r_oscillator = Relaxation(AlefeldPotraShi, x-> norm(x))

#sol_oscillator = solve(prob_oscillator, Tsit5_for_relaxation(); modif = r_oscillator)
#sol_exact = [prob_oscillator.f.analytic(prob_oscillator.u0, prob_oscillator.p, t) for t in sol_oscillator.t]
#plot(sol_oscillator)
#plot!(sol_oscillator.t, [sol_exact[i][1] for i ∈ 1:length(sol_oscillator.t)], label = "exact u[1]", lw = 4)
#plot!(sol_oscillator.t, [sol_exact[i][2] for i ∈ 1:length(sol_oscillator.t)], label = "exact u[2]", lw = 4)

sim = test_convergence(dts, prob_oscillator, Tsit5_for_relaxation())
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))
sim = test_convergence(dts, prob_oscillator, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim.𝒪est[:final]))
sim = test_convergence(dts, prob_oscillator, Tsit5_for_relaxation(); modif = r_oscillator)
println("order of convergence of new perform_step! with relaxation: "*string(sim.𝒪est[:final]))



########################################################
# TEST  2 : Non Linear Oscillator
printstyled("Non linear Harmonic Oscillator\n"; bold = true)

f_nloscillator = (u, p, t) -> [-u[2]/(u[1]^2 + u[2]^2),u[1]/(u[1]^2 + u[2]^2)]
prob_nloscillator = ODEProblem(
    ODEFunction(f_nloscillator; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))
r_nloscillator = Relaxation(AlefeldPotraShi, x-> norm(x))

sol_nloscillator = solve(prob_nloscillator, Tsit5_for_relaxation())
sol_exact = [prob_oscillator.f.analytic(prob_nloscillator.u0, prob_nloscillator.p, t) for t in sol_nloscillator.t]

sim = test_convergence(dts, prob_nloscillator, Tsit5_for_relaxation())
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))
sim = test_convergence(dts, prob_nloscillator, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim.𝒪est[:final]))
sim = test_convergence(dts, prob_nloscillator, Tsit5_for_relaxation(); modif = r_nloscillator)
println("order of convergence of new perform_step! with relaxation: "*string(sim.𝒪est[:final]))


########################################################
# TEST  3 : Non Linear Pendulum
printstyled("Non linear Pendulum\n"; bold = true)

f_nlpendulum = (u, p, t) -> [-sin(u[2]), u[1]]
prob_nlpendulum = ODEProblem(
    f_nlpendulum,
    [1.0, 0.0],
    (0.0, 1.0))
r_nlpendulum = Relaxation(AlefeldPotraShi, x-> x[1]^2/2 -  cos(x[2]))

#sol_nlpendulum = solve(prob_nlpendulum, Tsit5_for_relaxation(); modif = f_nlpendulum)
#sol_ref = solve(prob_nlpendulum, Vern9())

#=
sim = analyticless_test_convergence(dts, prob_nlpendulum, Tsit5_for_relaxation())
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))
sim = analyticless_test_convergence(dts, prob_nlpendulum, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim.𝒪est[:final]))
sim = analyticless_test_convergence(dts, prob_nlpendulum, Tsit5_for_relaxation(); modif = r_nlpendulum)
println("order of convergence of new perform_step! with relaxation: "*string(sim.𝒪est[:final]))
=#

=#