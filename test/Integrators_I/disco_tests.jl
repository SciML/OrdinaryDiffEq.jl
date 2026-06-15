using OrdinaryDiffEqCore, DiffEqDevTools, Test, LinearAlgebra
using OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqFIRK

predictive_disco_controller(alg) = OrdinaryDiffEqCore.PredictiveController(alg; discontinuity_detection = true)
PI_disco_controller(alg) = OrdinaryDiffEqCore.PIController(alg; discontinuity_detection = true)

function default_affect!(integrator)
    return nothing
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

sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6, controller = predictive_disco_controller(RadauIIA5()))
#  310.958 μs (8429 allocations: 263.88 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  392.667 μs (10628 allocations: 331.61 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  443.834 μs (15594 allocations: 547.55 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#  620.000 μs (21828 allocations: 765.69 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  60.250 μs (7576 allocations: 240.67 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  58.333 μs (7570 allocations: 240.58 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 2: TWO DISCONTINUITIES
#two discontinuity functions
function f(u, p, t)
    return if u[1] < 1
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

condition2(u, t, integrator) = u[1] - 2
cb2 = ContinuousCallback(condition2, default_affect!; maybe_discontinuity = true)
cb = CallbackSet(cb1, cb2)

sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6, controller = predictive_disco_controller(RadauIIA5()))
#  1.067 ms (28126 allocations: 875.64 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  1.311 ms (34747 allocations: 1.05 MiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  1.128 ms (38426 allocations: 1.31 MiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#  1.464 ms (52903 allocations: 1.80 MiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  271.125 μs (35982 allocations: 1.11 MiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  279.166 μs (39645 allocations: 1.22 MiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 3: EXPONENTIAL DISCONTINUITY
# multiple exponential regions with sharp transitions
function f_multi_exp!(du, u, p, t)
    return if u[1] < 0.3
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

cond_multi_2(u, t, integrator) = u[1] - 0.8
cb_multi_2 = ContinuousCallback(cond_multi_2, default_affect!; maybe_discontinuity = true)

cb_multi = CallbackSet(cb_multi_1, cb_multi_2)

sol_disco_radau = solve(prob_multi, RadauIIA5(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(RadauIIA5()))
#  137.542 μs (1470 allocations: 61.42 KiB)
sol_no_disco_radau = solve(prob_multi, RadauIIA5(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  127.209 μs (1073 allocations: 48.44 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Rodas5P()))
#  229.167 μs (1593 allocations: 64.94 KiB)
sol_no_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  239.750 μs (969 allocations: 45.44 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  96.917 μs (1469 allocations: 56.32 KiB)
sol_no_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol = 1e-7, abstol = 1e-9)
#  90.250 μs (1133 allocations: 46.84 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 4: STIFF MULTI-COMPONENT DISCONTINUITY
function f_stiff_disc!(du, u, p, t)
    λ = p[1]
    if u[1] < 0.5
        du[1] = -λ * u[1] + λ * exp(-t)
        du[2] = -λ * (u[2] - u[1])
        du[3] = -(u[3] - u[1]^2)
    else
        du[1] = u[1] - u[2]
        du[2] = -λ * (u[2] - u[1])
        du[3] = -(u[3] - u[1])
    end
end

u0_stiff = [0.1, 0.1, 0.01]
tspan_stiff = (0.0, 2.0)
prob_stiff = ODEProblem(f_stiff_disc!, u0_stiff, tspan_stiff, [500.0])

cond_stiff(u, t, integrator) = u[1] - 0.5
cb_stiff = ContinuousCallback(cond_stiff, default_affect!; maybe_discontinuity = true)

sol_disco_radau = solve(prob_stiff, RadauIIA5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11, controller = predictive_disco_controller(RadauIIA5()))
#  202.750 μs (2250 allocations: 101.73 KiB)
sol_no_disco_radau = solve(prob_stiff, RadauIIA5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  198.958 μs (2054 allocations: 95.58 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11, controller = PI_disco_controller(Rodas5P()))
#  379.875 μs (2053 allocations: 85.81 KiB)
sol_no_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  335.709 μs (1765 allocations: 79.78 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

#TEST 5: DISCONTINUOUS DAE
# discontinuous DAE with mass matrix
# System: M * du/dt = f(u, p, t)
# du[1]/dt = u[2] - u[1]
# 0 = u[1] + u[2] - 1 (algebraic constraint)
function f_dae_disc!(du, u, p, t)
    return if u[1] < 0.5
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

f_dae_func = ODEFunction(f_dae_disc!; mass_matrix = M_dae)
prob_dae = ODEProblem(f_dae_func, u0_dae, tspan_dae)

cond_dae(u, t, integrator) = u[1] - 0.5
cb_dae = ContinuousCallback(cond_dae, default_affect!; maybe_discontinuity = true)

sol_disco_rosenbrock = solve(prob_dae, Rodas5P(); callback=cb_dae, reltol=1e-8, abstol=1e-10, controller = PI_disco_controller(Rodas5P()))
#  117.500 μs (810 allocations: 35.44 KiB)
sol_no_disco_rosenbrock = solve(prob_dae, Rodas5P(); callback=cb_dae, reltol=1e-8, abstol=1e-10)
#  123.333 μs (625 allocations: 31.69 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

#TEST 6: VECTOR CALLBACK
function f_vec_disc!(du, u, p, t)
    du[1] = u[1] < 1.5 ? 2.0 * u[1] : -u[1] + 5.25
    du[2] = u[2] < 0.5 ? 3.0 * u[2] : -2.0 * u[2] + 2.5
    return du[3] = 4.0 * u[1] - u[2]
end

u0_vec = [0.1, 0.05, -0.1]
tspan_vec = (0.0, 3.0)
prob_vec = ODEProblem(f_vec_disc!, u0_vec, tspan_vec)

function condition_vec!(out, u, t, integrator)
    out[1] = u[1] - 1.5
    return out[2] = u[2] - 0.5
end
default_affect_vec!(integrator, idx) = nothing
cb_vec = VectorContinuousCallback(condition_vec!, default_affect_vec!, 2; maybe_discontinuity = true)

sol_disco_rosenbrock = solve(prob_vec, Rodas5P(); callback = cb_vec, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Rodas5P()))
#  168.250 μs (1109 allocations: 51.56 KiB)
sol_no_disco_rosenbrock = solve(prob_vec, Rodas5P(); callback = cb_vec, reltol=1e-7, abstol=1e-9)
#  204.292 μs (937 allocations: 49.41 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  98.542 μs (1705 allocations: 69.88 KiB)
sol_no_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol=1e-7, abstol=1e-9)
#  80.667 μs (1354 allocations: 62.98 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 7
function f!(du, u, p, t)
    x1, x2 = u
    du[1] = x2
    return if x2 < 0.0
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

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-8, abstol = 1e-10, controller = PI_disco_controller(Rodas5P()))
#  152.583 μs (1074 allocations: 48.30 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#  181.875 μs (877 allocations: 45.72 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10, controller = PI_disco_controller(Tsit5()))
#  62.791 μs (1284 allocations: 55.26 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#  54.708 μs (1155 allocations: 54.55 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject
