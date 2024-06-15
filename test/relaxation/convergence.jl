using OrdinaryDiffEq, DiffEqDevTools,  Test

import ODEProblemLibrary: prob_ode_linear,
                          prob_ode_2Dlinear

include("relaxation.jl")

#########################################################################
# Check that the code with the new structure without modification asked 
# by the user gives the same result
# Comparaison on Tsit5

probnum = prob_ode_linear
prob = prob_ode_2Dlinear


testTol = 0.2

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


#########################################################################
##                            With Relaxation 

dts = (1 / 2) .^ (6:-1:4)

########################################################
# TEST  1 : Harmonic Oscillator
printstyled("Harmonic Oscillator\n"; bold = true)

f_oscillator = (u, p, t) -> [-u[2],u[1]]
prob_oscillator = ODEProblem(
    ODEFunction(f_oscillator; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))
r_oscillator = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> norm(x)))

#sol_oscillator = solve(prob_oscillator, Tsit5_for_relaxation(); modif = r_oscillator)
#sol_exact = [prob_oscillator.f.analytic(prob_oscillator.u0, prob_oscillator.p, t) for t in sol_oscillator.t]
#plot(sol_oscillator)
#plot!(sol_oscillator.t, [sol_exact[i][1] for i ∈ 1:length(sol_oscillator.t)], label = "exact u[1]", lw = 4)
#plot!(sol_oscillator.t, [sol_exact[i][2] for i ∈ 1:length(sol_oscillator.t)], label = "exact u[2]", lw = 4)

sim_oscillator_old = test_convergence(dts, prob_oscillator, Tsit5())
println("order of convergence of older perform_step! : "*string(sim_oscillator_old.𝒪est[:final]))
sim_oscillator_new = test_convergence(dts, prob_oscillator, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim_oscillator_new.𝒪est[:final]))
sim_oscillator_relax = test_convergence(dts, prob_oscillator, Tsit5_for_relaxation(); modif = r_oscillator)
println("order of convergence of new perform_step! with relaxation: "*string(sim_oscillator_relax.𝒪est[:final]))

########################################################
# TEST  2 : Non Linear Oscillator
printstyled("Non linear Harmonic Oscillator\n"; bold = true)

f_nloscillator = (u, p, t) -> [-u[2]/(u[1]^2 + u[2]^2),u[1]/(u[1]^2 + u[2]^2)]
prob_nloscillator = ODEProblem(
    ODEFunction(f_nloscillator; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))
r_nloscillator = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> norm(x)))

#sol_nloscillator = solve(prob_nloscillator, Tsit5_for_relaxation(); modif = r_nloscillator)
#sol_exact = [prob_oscillator.f.analytic(prob_nloscillator.u0, prob_nloscillator.p, t) for t in sol_nloscillator.t]

sim = test_convergence(dts, prob_nloscillator, Tsit5())
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
r_nlpendulum = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> x[1]^2/2 -  cos(x[2])))

#sol_nlpendulum = solve(prob_nlpendulum, Tsit5_for_relaxation(); modif = r_nlpendulum)
#sol_ref = solve(prob_nlpendulum, Vern9())

test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)

sim = analyticless_test_convergence(dts, prob_nlpendulum, Tsit5(), test_setup)
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))
sim = analyticless_test_convergence(dts, prob_nlpendulum, Tsit5_for_relaxation(), test_setup)
println("order of convergence of new perform_step! without relaxation: "*string(sim.𝒪est[:final]))
sim = analyticless_test_convergence(dts, prob_nlpendulum, Tsit5_for_relaxation(), test_setup; modif = r_nlpendulum)
println("order of convergence of new perform_step! with relaxation: "*string(sim.𝒪est[:final]))


@benchmark solve(prob_nloscillator, Tsit5_for_relaxation())

############################################################################
# TEST  4 : Time dependent harmonic oscillator with bounded angular velocity

printstyled("Time dependent harmonic oscillator\n"; bold = true)

f_td_oscillator = (u, p, t) -> [-(1+0.5 * sin(t))*u[2], (1+0.5 * sin(t))*u[1]]
prob_td_oscillator = ODEProblem(
    ODEFunction(f_td_oscillator; 
                analytic = (u0, p, t)->[cos(0.5)*cos(t-0.5*cos(t))-sin(0.5)*sin(t-0.5*cos(t)), 
                                        sin(0.5)*cos(t-0.5*cos(t))+cos(0.5)*sin(t-0.5*cos(t))]),
    [1.0, 0.0],
    (0.0, 1.0))
r_td_oscillator = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> norm(x)))

sim_td_oscillator_old = test_convergence(dts, prob_td_oscillator, Tsit5())
println("order of convergence of older perform_step! : "*string(sim_td_oscillator_old.𝒪est[:final]))
sim_td_oscillator_new = test_convergence(dts, prob_td_oscillator, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim_td_oscillator_new.𝒪est[:final]))
sim_td_oscillator_relax = test_convergence(dts, prob_td_oscillator, Tsit5_for_relaxation(); modif = r_oscillator)
println("order of convergence of new perform_step! with relaxation: "*string(sim_td_oscillator_relax.𝒪est[:final]))

############################################################################
# TEST  5 : Conserved exponential entropy

printstyled("Conserved exponential entropy\n"; bold = true)

f_cee = (u, p, t) -> [-exp(u[2]), exp(u[1])]
prob_cee  = ODEProblem(
    ODEFunction(f_cee ; 
                analytic = (u0, p, t)->[log(exp(1) + exp(0.5)) - log(exp(0.5) + exp((exp(0.5)+exp(1))*t)), 
                                        log(exp((exp(0.5)+exp(1))*t)*(exp(0.5)+exp(1)))/(exp(0.5) + exp((exp(0.5)+exp(1))*t))]),
    [1.0, 0.5],
    (0.0, 1.0))
r_cee  = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> exp(x[1]) + exp(x[2])))

sim_cee_old = test_convergence(dts, prob_cee, Tsit5())
println("order of convergence of older perform_step! : "*string(sim_cee_old.𝒪est[:final]))
sim_cee_new = test_convergence(dts, prob_cee, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim_cee_new.𝒪est[:final]))
sim_cee_relax = test_convergence(dts, prob_cee, Tsit5_for_relaxation(); modif = r_oscillator)
println("order of convergence of new perform_step! with relaxation: "*string(sim_cee_relax.𝒪est[:final]))