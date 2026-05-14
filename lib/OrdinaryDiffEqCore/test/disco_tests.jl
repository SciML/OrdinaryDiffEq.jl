using OrdinaryDiffEqCore, DiffEqDevTools, Test, LinearAlgebra
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqFIRK
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error)) 
using BenchmarkTools

predictive_disco_controller(alg) = OrdinaryDiffEqCore.PredictiveController(alg; discontinuity_detection = true)
PI_disco_controller(alg) = OrdinaryDiffEqCore.PIController(alg; discontinuity_detection = true)

function default_affect!(integrator)
    nothing
end

#TEST 1: SIMPLE DISCONTINUITY
#test example discontinuous at u = 1
f(u, p, t) = u[1] < 1 ? [2u[1]] : [-3u[1] + 5]
u0 = [0.1]           
tspan = (0.0, 1.5)
prob = ODEProblem(f, u0, tspan)

#define callback
condition(u, t, integrator) = u[1] - 1
cb = ContinuousCallback(condition, default_affect!; maybe_discontinuity = true)
cb2 = ContinuousCallback(condition, default_affect!; maybe_discontinuity = false)

sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6, controller = predictive_disco_controller(RadauIIA5()))
#  251.334 μs (7692 allocations: 240.09 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb2, reltol = 1e-6)
#  339.167 μs (10308 allocations: 320.84 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  439.417 μs (16844 allocations: 589.44 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1e-6)
#  580.125 μs (21804 allocations: 763.61 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  49.917 μs (7180 allocations: 227.70 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1e-6)
#  49.542 μs (7158 allocations: 227.05 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 2: TWO DISCONTINUITIES
#two discontinuity functions
function f(u, p, t)
    if u[1] < 1
        [2u[1]]                
    elseif u[1] < 2
        [u[1] + 0.2]           
    else
        [-4u[1] + 12]           
    end
end

u0 = [0.1]
tspan = (0.0, 2.5)
prob = ODEProblem(f, u0, tspan)

#define callbacks
condition1(u, t, integrator) = u[1] - 1
cb1 = ContinuousCallback(condition1, default_affect!; maybe_discontinuity = true)
cb1f = ContinuousCallback(condition1, default_affect!; maybe_discontinuity = false)

condition2(u, t, integrator) = u[1] - 2
cb2 = ContinuousCallback(condition2, default_affect!; maybe_discontinuity = true)
cb2f = ContinuousCallback(condition2, default_affect!; maybe_discontinuity = false)
cb = CallbackSet(cb1, cb2)
cb2 = CallbackSet(cb1f, cb2f)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  1.279 ms (45803 allocations: 1.56 MiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1e-6)
#  1.412 ms (52871 allocations: 1.80 MiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  274.709 μs (34714 allocations: 1.06 MiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1e-6)
#  259.625 μs (38745 allocations: 1.19 MiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 3: EXPONENTIAL DISCONTINUITY
# multiple exponential regions with sharp transitions
function f_multi_exp!(du, u, p, t)
    if u[1] < 0.3
        du[1] = 3 * exp(3 * u[1])    
    elseif u[1] < 0.8
        du[1] = exp(u[1])    
    else
        du[1] = u[1]  
    end
end

u0_multi = [0.05]
tspan_multi = (0.0, 1.5)
prob_multi = ODEProblem(f_multi_exp!, u0_multi, tspan_multi)

#define callbacks
cond_multi_1(u, t, integrator) = u[1] - 0.3
cb_multi_1 = ContinuousCallback(cond_multi_1, default_affect!; maybe_discontinuity = true)
cb_multi_1f = ContinuousCallback(cond_multi_1, default_affect!; maybe_discontinuity = false)

cond_multi_2(u, t, integrator) = u[1] - 0.8
cb_multi_2 = ContinuousCallback(cond_multi_2, default_affect!; maybe_discontinuity = true)
cb_multi_2f = ContinuousCallback(cond_multi_2, default_affect!; maybe_discontinuity = false)

cb_multi = CallbackSet(cb_multi_1, cb_multi_2)
cb_multi2 = CallbackSet(cb_multi_1f, cb_multi_2f)

#disco solve
sol_disco_radau = solve(prob_multi, RadauIIA5(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = predictive_disco_controller(RadauIIA5()))
#  164.000 μs (2046 allocations: 76.33 KiB)
sol_no_disco_radau = solve(prob_multi, RadauIIA5(); callback=cb_multi2, reltol = 1e-7, abstol = 1e-9)
#  127.541 μs (1044 allocations: 46.11 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Rodas5P()))
#  308.584 μs (2913 allocations: 103.64 KiB)
sol_no_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi2, reltol=1e-7, abstol=1e-9)
#  236.833 μs (940 allocations: 43.33 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  120.167 μs (2240 allocations: 83.23 KiB)
sol_no_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi2, reltol = 1e-7, abstol = 1e-9)
#  88.833 μs (1105 allocations: 45.27 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

@profview for i in 1:1000 solve(prob_multi, Tsit5(); callback = cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5())) end

#TEST 4: STIFF DISCONTINUITY
# very stiff discontinuous system
function f_stiff_disc!(du, u, p, t)
    λ = p[1] # stiffness parameter
    if u[1] < 0.5
        du[1] = -λ * u[1] + λ * exp(-t) # stiff decay with forcing
    else
        du[1] = u[1]
    end
end

u0_stiff = [0.1]
tspan_stiff = (0.0, 3.0)
prob_stiff = ODEProblem(f_stiff_disc!, u0_stiff, tspan_stiff, [100.0])

#define callback
cond_stiff(u, t, integrator) = u[1] - 0.5
cb_stiff = ContinuousCallback(cond_stiff, default_affect!; maybe_discontinuity = true)
cb_stiff_f = ContinuousCallback(cond_stiff, default_affect!; maybe_discontinuity = false)

#disco solve
sol_disco_radau = solve(prob_stiff, RadauIIA5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11, controller = predictive_disco_controller(RadauIIA5()))
#  123.708 μs (1786 allocations: 69.94 KiB)
sol_no_disco_radau = solve(prob_stiff, RadauIIA5(); callback=cb_stiff_f, reltol = 1e-9, abstol = 1e-11)
#  125.416 μs (1476 allocations: 61.14 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_tsit5 = solve(prob_stiff, Tsit5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11, controller = PI_disco_controller(Tsit5()))
#  85.208 μs (2006 allocations: 75.09 KiB)
sol_no_disco_tsit5 = solve(prob_stiff, Tsit5(); callback=cb_stiff_f, reltol = 1e-9, abstol = 1e-11)
#  80.833 μs (1783 allocations: 68.51 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 5: DISCONTINUOUS DAE
# discontinuous DAE with mass matrix
# System: M * du/dt = f(u, p, t)
# du[1]/dt = u[2] - u[1]
# 0 = u[1] + u[2] - 1 (algebraic constraint)
function f_dae_disc!(du, u, p, t)
    if u[1] < 0.5
        du[1] = 2 * u[2] - u[1]
        du[2] = u[1] + u[2] - 1 # algebraic constraint 
    else
        du[1] = -u[1] + u[2]
        du[2] = u[1] + u[2] - 1  
    end
end

u0_dae = [0.2, 0.8]
tspan_dae = (0.0, 2.0)

M_dae = [1.0 0.0; 0.0 0.0]

f_dae_func = ODEFunction(f_dae_disc!; mass_matrix=M_dae)
prob_dae = ODEProblem(f_dae_func, u0_dae, tspan_dae)

cond_dae(u, t, integrator) = u[1] - 0.5
cb_dae = ContinuousCallback(cond_dae, default_affect!; maybe_discontinuity = true)
cb_daef = ContinuousCallback(cond_dae, default_affect!; maybe_discontinuity = false)

radau_disco = solve(prob_dae, RadauIIA5(); callback=cb_dae, reltol=1e-8, abstol=1e-10, controller = predictive_disco_controller(RadauIIA5()))
#  66.791 μs (827 allocations: 36.17 KiB)
radau_no_disco = solve(prob_dae, RadauIIA5(); callback=cb_daef, reltol=1e-8, abstol=1e-10)
#  61.166 μs (639 allocations: 30.59 KiB)
@test radau_disco.retcode == ReturnCode.Success
@test radau_disco.stats.nreject <= radau_no_disco.stats.nreject

#TEST 6: VECTOR CALLBACK
function f!(du, u, p, t)
    du[1] = -u[1]
    du[2] =  0.2*u[1] - 0.1*u[2]
    du[3] = 0.1*u[2] - 0.1*u[3]
end

u0    = [3.0, 0.0, -2.0]
tspan = (0.0, 10.0)
prob  = ODEProblem(f!, u0, tspan)

# u[1] == 2.0 and u[1] == 1.0
function condition!(out, u, t, integrator)
    out[1] = u[1] - 2.0
    out[2] = u[1] - 1.0
end

# Discontinuous update to the state when an event fires
function affect!(integrator, idx)
    if idx == 1
        # when u[1] crosses 2, kick u[2] up (jump discontinuity)
        integrator.u[2] += 5.0
    elseif idx == 2
        # when u[1] crosses 1, reset u[2]
        integrator.u[2] = 0.0
    end
end

cb = VectorContinuousCallback(condition!, affect!, 2; maybe_discontinuity = true)
cb2 = VectorContinuousCallback(condition!, affect!, 2; maybe_discontinuity = false)  

sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, controller = predictive_disco_controller(RadauIIA5()))
#  29.500 μs (476 allocations: 23.75 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb2)
#  29.292 μs (465 allocations: 23.09 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, controller = PI_disco_controller(Rodas5P()))
#  44.959 μs (412 allocations: 23.12 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2)
#  44.291 μs (401 allocations: 22.56 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, controller = PI_disco_controller(Tsit5()))
#  23.208 μs (508 allocations: 23.54 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2)
#  23.375 μs (497 allocations: 22.91 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 7
function f!(du, u, p, t)
    x1, x2 = u
    du[1] = x2
    if x2 < 0.0
        du[2] = -x1 + 1.0
    else
        du[2] = -x1 - 1.0
    end
end

u = [1.5, 0.8]
tspan = (0.0, 2.0)
prob = ODEProblem(f!, u, tspan)

cond(u, t, integrator) = u[2]
cb = ContinuousCallback(cond, default_affect!; maybe_discontinuity = true)
cb2 = ContinuousCallback(cond, default_affect!; maybe_discontinuity = false)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-8, abstol = 1e-10, controller = PI_disco_controller(Rodas5P()))
#  245.875 μs (2496 allocations: 84.00 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1e-8, abstol = 1e-10)
#  174.291 μs (845 allocations: 43.62 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10, controller = PI_disco_controller(Tsit5()))
#  76.583 μs (2053 allocations: 75.30 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1e-8, abstol = 1e-10)
#  53.250 μs (1117 allocations: 52.55 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject