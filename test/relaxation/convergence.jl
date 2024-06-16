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

########################################################
# TEST  1 : Harmonic Oscillator
include("harmonic_oscillator.jl")

########################################################
# TEST  2 : Non Linear Oscillator
include("non_linear_oscillator.jl")

########################################################
# TEST  3 : Non Linear Pendulum
include("non_linear_pendulum.jl")

############################################################################
# TEST  4 : Time dependent harmonic oscillator with bounded angular velocity
include("time_dependant_harmonic_oscillator.jl")

############################################################################
# TEST  5 : Conserved exponential entropy
include("conserved_exponential_entropy.jl")