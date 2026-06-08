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
#  272.333 μs (7961 allocations: 249.34 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  369.000 μs (10315 allocations: 321.41 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  419.750 μs (15586 allocations: 546.83 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#  587.042 μs (21827 allocations: 764.73 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  54.916 μs (7185 allocations: 228.16 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  53.083 μs (7165 allocations: 227.38 KiB)
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
#  917.292 μs (25563 allocations: 795.66 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  1.325 ms (36765 allocations: 1.11 MiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Rodas5P()))
#  1.371 ms (48673 allocations: 1.66 MiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-6)
#  1.465 ms (52921 allocations: 1.80 MiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6, controller = PI_disco_controller(Tsit5()))
#  272.292 μs (34890 allocations: 1.07 MiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-6)
#  285.625 μs (38807 allocations: 1.19 MiB)
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

sol_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Rodas5P()))
#  210.333 μs (1403 allocations: 58.92 KiB)
sol_no_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#  238.250 μs (1024 allocations: 45.09 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  98.084 μs (1488 allocations: 58.20 KiB)
sol_no_disco_tsit5 = solve(prob_multi, Tsit5(); callback=cb_multi, reltol = 1e-7, abstol = 1e-9)
#  89.958 μs (1184 allocations: 47.15 KiB)
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
#  208.959 μs (2257 allocations: 101.69 KiB)
sol_no_disco_radau = solve(prob_stiff, RadauIIA5(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  191.500 μs (2009 allocations: 93.05 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11, controller = PI_disco_controller(Rodas5P()))
#  340.958 μs (2060 allocations: 105.09 KiB)
sol_no_disco_rosenbrock = solve(prob_stiff, Rodas5P(); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#  325.375 μs (1793 allocations: 79.44 KiB)
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

radau_disco = solve(prob_dae, RadauIIA5(); callback=cb_dae, reltol=1e-8, abstol=1e-10, controller = predictive_disco_controller(RadauIIA5()))
#  62.083 μs (743 allocations: 33.97 KiB)
radau_no_disco = solve(prob_dae, RadauIIA5(); callback=cb_dae, reltol=1e-8, abstol=1e-10)
#  58.458 μs (641 allocations: 30.88 KiB)
@test radau_disco.retcode == ReturnCode.Success
@test radau_disco.stats.nreject <= radau_no_disco.stats.nreject

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
#  197.417 μs (1289 allocations: 57.34 KiB)
sol_no_disco_rosenbrock = solve(prob_vec, Rodas5P(); callback = cb_vec, reltol=1e-7, abstol=1e-9)
#  192.542 μs (936 allocations: 47.69 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol=1e-7, abstol=1e-9, controller = PI_disco_controller(Tsit5()))
#  87.583 μs (1537 allocations: 63.20 KiB)
sol_no_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol=1e-7, abstol=1e-9)
#  76.042 μs (1295 allocations: 59.73 KiB)
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
#  161.167 μs (1154 allocations: 51.39 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#  177.250 μs (910 allocations: 45.09 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10, controller = PI_disco_controller(Tsit5()))
#  64.958 μs (1357 allocations: 59.46 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1e-8, abstol = 1e-10)
#  52.375 μs (1173 allocations: 54.07 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject
