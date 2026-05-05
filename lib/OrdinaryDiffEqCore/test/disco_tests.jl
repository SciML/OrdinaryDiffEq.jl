using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqLowOrderRK
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error)) 
using BenchmarkTools

#TEST 1: SIMPLE DISCONTINUITY
#test example discontinuous at u = 1
f(u, p, t) = u[1] < 1 ? [2u[1]] : [-3u[1] + 5]
u0 = [0.1]           
tspan = (0.0, 1.5)
prob = ODEProblem(f, u0, tspan)

#define callback
condition(u, t, integrator) = u[1] - 1
function affect!(integrator)
    #println("fired callback at t=$(integrator.t), u=$(integrator.u[1])")
    integrator.u[1] += 10
end
cb = ContinuousCallback(condition, affect!; is_discontinuity = true)
cb2 = ContinuousCallback(condition, affect!; is_discontinuity = false)

sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  294.458 μs (8082 allocations: 256.59 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb2, reltol = 1e-6)
#  356.708 μs (10024 allocations: 312.08 KiB)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#   474.375 μs (16801 allocations: 592.50 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1e-6)
#  509.083 μs (18240 allocations: 639.33 KiB)

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  59.542 μs (7248 allocations: 233.67 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1e-6)
#  46.500 μs (7129 allocations: 226.22 KiB)

sol_disco_BS5 = solve(prob, BS5(); callback = cb, reltol = 1e-6)
sol_no_disco_BS5 = solve(prob, BS5(); callback = cb2, reltol = 1e-6)

#TEST 2: TWO DISCONTINUITIES
#two discontinuity functions
function f(u, p, t)
    if u[1] < 1
        [2u[1]]                 # region 1: grows to hit u = 1
    elseif u[1] < 2
        [u[1] + 0.2]            # region 2: continues increasing to hit u = 2
    else
        [-4u[1] + 12]           
    end
end

u0 = [0.1]
tspan = (0.0, 2.5)
prob = ODEProblem(f, u0, tspan)

#define callbacks
condition1(u, t, integrator) = u[1] - 1
function affect1!(integrator)
    #println("Callback 1 fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb1 = ContinuousCallback(condition1, affect1!; is_discontinuity = true)
cb1f = ContinuousCallback(condition1, affect1!; is_discontinuity = false)

condition2(u, t, integrator) = u[1] - 2
function affect2!(integrator)
    #println("Callback 2 fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb2 = ContinuousCallback(condition2, affect2!; is_discontinuity = true)
cb2f = ContinuousCallback(condition2, affect2!; is_discontinuity = false)
cb = CallbackSet(cb1, cb2)
cb2 = CallbackSet(cb1f, cb2f)

#disco solve
sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  1.503 ms (41672 allocations: 1.27 MiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb2, reltol = 1e-6)
#  1.306 ms (37092 allocations: 1.13 MiB)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#   1.164 ms (44318 allocations: 1.52 MiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1e-6)
#   1.306 ms (51713 allocations: 1.76 MiB)

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#   279.792 μs (34573 allocations: 1.07 MiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1e-6)
#  266.167 μs (39024 allocations: 1.21 MiB)


#TEST 3: EXPONENTIAL DISCONTINUITY
# multiple exponential regions with sharp transitions
function f_multi_exp!(du, u, p, t)
    if u[1] < 0.3
        du[1] = 3 * exp(3 * u[1])    # very steep exponential
    elseif u[1] < 0.8
        du[1] = exp(u[1])      # slower exponential
    else
        du[1] = u[1]  # linear
    end
end

u0_multi = [0.05]
tspan_multi = (0.0, 1.5)
prob_multi = ODEProblem(f_multi_exp!, u0_multi, tspan_multi)

#define callbacks
cond_multi_1(u, t, integrator) = u[1] - 0.3
function affect_multi_1!(integrator)
    #println("Multi-exponential discontinuity 1 callback fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb_multi_1 = ContinuousCallback(cond_multi_1, affect_multi_1!; is_discontinuity = true)
cb_multi_1f = ContinuousCallback(cond_multi_1, affect_multi_1!; is_discontinuity = false)

cond_multi_2(u, t, integrator) = u[1] - 0.8
function affect_multi_2!(integrator)
    #println("Multi-exponential discontinuity 2 callback fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb_multi_2 = ContinuousCallback(cond_multi_2, affect_multi_2!; is_discontinuity = true)
cb_multi_2f = ContinuousCallback(cond_multi_2, affect_multi_2!; is_discontinuity = false)
cb_multi = CallbackSet(cb_multi_1, cb_multi_2)
cb_multi2 = CallbackSet(cb_multi_1f, cb_multi_2f)

#disco solve
sol_disco = solve(prob_multi, RadauIIA5(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  175.625 μs (1871 allocations: 81.55 KiB)
sol_no_disco = solve(prob_multi, RadauIIA5(); callback=cb_multi2, reltol = 1e-7, abstol = 1e-9)
#  142.875 μs (1244 allocations: 59.17 KiB)

sol_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#   295.834 μs (2216 allocations: 90.70 KiB)
sol_no_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi2, reltol=1e-7, abstol=1e-9)
#  253.709 μs (1380 allocations: 74.28 KiB)

sol_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  127.375 μs (1953 allocations: 87.49 KiB)
sol_no_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi2, reltol = 1e-7, abstol = 1e-9)
#  95.250 μs (1499 allocations: 73.62 KiB)

sol_disco_BS3 = solve(prob_multi, BS3(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
sol_no_disco_BS3 = solve(prob_multi, BS3(); callback=cb_multi2, reltol=1e-7, abstol=1e-9)

@profview for i in 1:1000 
    solve(prob_multi, RadauIIA5(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
end

#TEST 4: STIFF DISCONTINUITY
# very stiff discontinuous system
function f_stiff_disc!(du, u, p, t)
    λ = p[1]  # stiffness parameter
    if u[1] < 0.5
        du[1] = -λ * u[1] + λ * exp(-t)  # stiff decay with forcing
    else
        du[1] = u[1]
    end
end

u0_stiff = [0.1]
tspan_stiff = (0.0, 3.0)
prob_stiff = ODEProblem(f_stiff_disc!, u0_stiff, tspan_stiff, [100.0])

#define callback
cond_stiff(u, t, integrator) = u[1] - 0.5
function affect_stiff!(integrator)
    #println("Stiff discontinuity callback fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb_stiff = ContinuousCallback(cond_stiff, affect_stiff!; is_discontinuity = true)
cb_stiff_f = ContinuousCallback(cond_stiff, affect_stiff!; is_discontinuity = false)

#disco solve
sol_disco = solve(prob_stiff, RadauIIA5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  149.167 μs (1819 allocations: 75.19 KiB)
sol_no_disco = solve(prob_stiff, RadauIIA5(); callback=cb_stiff_f, reltol = 1e-9, abstol = 1e-11)
#  138.125 μs (1565 allocations: 64.09 KiB)

sol_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#   204.833 μs (1517 allocations: 59.33 KiB)
sol_no_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff_f, reltol=1e-9, abstol=1e-11)
#   156.500 μs (1047 allocations: 44.59 KiB)

sol_disco_tsit5 = solve(prob_stiff, Tsit5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  93.833 μs (2040 allocations: 80.59 KiB)
sol_no_disco_tsit5 = solve(prob_stiff, Tsit5(); callback=cb_stiff_f, reltol = 1e-9, abstol = 1e-11)
#  82.750 μs (1898 allocations: 72.12 KiB)

sol_disco_BS3 = solve(prob_stiff, BS3(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#   1.121 ms (12460 allocations: 595.30 KiB)
sol_no_disco_BS3 = solve(prob_stiff, BS3(); callback=cb_stiff_f, reltol=1e-9, abstol=1e-11)
#   1.102 ms (12229 allocations: 582.34 KiB)

#TEST 5: DISCONTINUOUS DAE
# discontinuous DAE with mass matrix
# System: M * du/dt = f(u, p, t)
# du[1]/dt = u[2] - u[1]
# 0 = u[1] + u[2] - 1 (algebraic constraint)
function f_dae_disc!(du, u, p, t)
    if u[1] < 0.5
        du[1] = 2 * u[2] - u[1]
        du[2] = u[1] + u[2] - 1  # algebraic constraint 
    else
        du[1] = -u[1] + u[2]
        du[2] = u[1] + u[2] - 1  # algebraic constraint 
    end
end

u0_dae = [0.2, 0.8]  # consistent with constraint u[1] + u[2] = 1
tspan_dae = (0.0, 2.0)

M_dae = [1.0 0.0; 0.0 0.0]

f_dae_func = ODEFunction(f_dae_disc!; mass_matrix=M_dae)
prob_dae = ODEProblem(f_dae_func, u0_dae, tspan_dae)

cond_dae(u, t, integrator) = u[1] - 0.5
function affect_dae!(integrator)
    #println("DAE discontinuity callback fired at t=$(integrator.t), u=$(integrator.u)")
end
cb_dae = ContinuousCallback(cond_dae, affect_dae!; is_discontinuity = true)
cb_daef = ContinuousCallback(cond_dae, affect_dae!; is_discontinuity = false)

radau_disco = solve(prob_dae, RadauIIA5(); callback=cb_dae, reltol=1e-8, abstol=1e-10)
#  88.542 μs (870 allocations: 41.86 KiB)
radau_no_disco = solve(prob_dae, RadauIIA5(); callback=cb_daef, reltol=1e-8, abstol=1e-10)
#  73.000 μs (673 allocations: 32.05 KiB)

sol_disco_rosenbrock = solve(prob_dae, Rodas5P(); callback=cb_dae, reltol=1e-8, abstol=1e-10)
#   312.167 μs (1200 allocations: 48.73 KiB)
sol_no_disco_rosenbrock = solve(prob_dae, Rodas5P(); callback=cb_daef, reltol=1e-8, abstol=1e-10)
#   256.792 μs (672 allocations: 32.56 KiB)

#TEST 6: VECTOR CALLBACK
function f!(du, u, p, t)
    du[1] = -u[1]
    du[2] =  0.2*u[1] - 0.1*u[2]
end

u0    = [3.0, 0.0]
tspan = (0.0, 10.0)
prob  = ODEProblem(f!, u0, tspan)

# Two event surfaces: u[1] == 2.0 and u[1] == 1.0
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

cb = VectorContinuousCallback(condition!, affect!, 2; is_discontinuity = true)
cb2 = VectorContinuousCallback(condition!, affect!, 2; is_discontinuity = false)  

sol_disco = solve(prob, RadauIIA5(); callback = cb)
#   49.125 μs (664 allocations: 32.89 KiB)
sol_no_disco = solve(prob, RadauIIA5(); callback = cb2)
#   37.375 μs (531 allocations: 25.23 KiB)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb)
#   57.333 μs (592 allocations: 31.23 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2)
#  44.250 μs (476 allocations: 23.73 KiB)

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb)
#   37.833 μs (673 allocations: 31.80 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2)
#  24.958 μs (557 allocations: 24.23 KiB)

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
affect!(integrator) = nothing

cb = ContinuousCallback(cond, affect!; is_discontinuity = true)
cb2 = ContinuousCallback(cond, affect!; is_discontinuity = false)

sol_disco = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-8, abstol = 1e-10)
sol_no_disco = solve(prob, RadauIIA5(); callback = cb2, reltol = 1e-8, abstol = 1e-10)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#   240.291 μs (1821 allocations: 71.56 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1e-8, abstol = 1e-10)
#   184.625 μs (1029 allocations: 49.23 KiB)

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#   79.791 μs (1678 allocations: 73.85 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1e-8, abstol = 1e-10)
#   55.958 μs (1259 allocations: 57.04 KiB)