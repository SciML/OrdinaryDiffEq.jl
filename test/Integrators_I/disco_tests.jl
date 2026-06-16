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
#  314.500 μs (8909 allocations: 278.73 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  376.125 μs (10628 allocations: 331.61 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  444.500 μs (15596 allocations: 547.62 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#  620.709 μs (21830 allocations: 765.77 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  55.625 μs (7576 allocations: 240.67 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  55.083 μs (7570 allocations: 240.58 KiB)
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
#  930.334 μs (25771 allocations: 802.09 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  1.256 ms (34747 allocations: 1.05 MiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  1.224 ms (42698 allocations: 1.46 MiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#  1.551 ms (52905 allocations: 1.80 MiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  267.167 μs (34737 allocations: 1.07 MiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  279.416 μs (39645 allocations: 1.22 MiB)
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
#  148.875 μs (1564 allocations: 64.56 KiB)
sol_no_disco_radau = solve(prob_multi, RadauIIA5(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  129.875 μs (1073 allocations: 48.44 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Rodas5P()))
#  209.000 μs (1474 allocations: 60.94 KiB)
sol_no_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  237.250 μs (969 allocations: 45.44 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  105.250 μs (1592 allocations: 60.34 KiB)
sol_no_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol = 1e-7, abstol = 1e-9)
#  88.000 μs (1133 allocations: 46.84 KiB)
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
#  207.500 μs (2285 allocations: 101.72 KiB)
sol_no_disco_radau = solve(prob_stiff, RadauIIA5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  196.458 μs (2054 allocations: 95.58 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11, controller = PI_disco_controller(Rodas5P()))
#  364.833 μs (2030 allocations: 85.34 KiB)
sol_no_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  332.750 μs (1765 allocations: 79.78 KiB)
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
#  127.083 μs (831 allocations: 35.53 KiB)
sol_no_disco_rosenbrock = solve(prob_dae, Rodas5P(); callback=cb_dae, reltol=1e-8, abstol=1e-10)
#  121.625 μs (625 allocations: 31.69 KiB)
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
#  191.250 μs (1212 allocations: 53.62 KiB)
sol_no_disco_rosenbrock = solve(prob_vec, Rodas5P(); callback = cb_vec, reltol=1e-7, abstol=1e-9)
#  202.000 μs (939 allocations: 49.47 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  94.459 μs (1626 allocations: 67.73 KiB)
sol_no_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol=1e-7, abstol=1e-9)
#  93.084 μs (1356 allocations: 63.04 KiB)
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
#  169.208 μs (1190 allocations: 51.78 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#  184.792 μs (877 allocations: 45.72 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10, controller = PI_disco_controller(Tsit5()))
#  59.625 μs (1248 allocations: 53.79 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#  54.041 μs (1155 allocations: 54.55 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject
