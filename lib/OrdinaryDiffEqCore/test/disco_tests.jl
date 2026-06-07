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
cb2 = ContinuousCallback(condition, default_affect!; maybe_discontinuity = false)

sol_disco_radau = solve(prob, RadauIIA5(); callback = cb, reltol = 1.0e-6, controller = predictive_disco_controller(RadauIIA5()))
#  253.208 μs (7657 allocations: 239.47 KiB)
sol_no_disco_radau = solve(prob, RadauIIA5(); callback = cb2, reltol = 1.0e-6)
#  347.917 μs (10308 allocations: 320.86 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1.0e-6, controller = PI_disco_controller(Rodas5P()))
#  456.292 μs (16771 allocations: 587.12 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1.0e-6)
#  589.666 μs (21804 allocations: 763.61 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1.0e-6, controller = PI_disco_controller(Tsit5()))
#  48.917 μs (7173 allocations: 227.61 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1.0e-6)
#  48.750 μs (7158 allocations: 227.05 KiB)
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
cb1f = ContinuousCallback(condition1, default_affect!; maybe_discontinuity = false)

condition2(u, t, integrator) = u[1] - 2
cb2 = ContinuousCallback(condition2, default_affect!; maybe_discontinuity = true)
cb2f = ContinuousCallback(condition2, default_affect!; maybe_discontinuity = false)
cb = CallbackSet(cb1, cb2)
cb2 = CallbackSet(cb1f, cb2f)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1.0e-6, controller = PI_disco_controller(Rodas5P()))
#  1.306 ms (45129 allocations: 1.54 MiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1.0e-6)
#  1.500 ms (52871 allocations: 1.80 MiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1.0e-6, controller = PI_disco_controller(Tsit5()))
#  261.375 μs (34373 allocations: 1.06 MiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1.0e-6)
#  266.125 μs (38745 allocations: 1.19 MiB)
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
cb_multi_1f = ContinuousCallback(cond_multi_1, default_affect!; maybe_discontinuity = false)

cond_multi_2(u, t, integrator) = u[1] - 0.8
cb_multi_2 = ContinuousCallback(cond_multi_2, default_affect!; maybe_discontinuity = true)
cb_multi_2f = ContinuousCallback(cond_multi_2, default_affect!; maybe_discontinuity = false)

cb_multi = CallbackSet(cb_multi_1, cb_multi_2)
cb_multi2 = CallbackSet(cb_multi_1f, cb_multi_2f)

#disco solve
sol_disco_radau = solve(prob_multi, RadauIIA5(); callback = cb_multi, reltol = 1.0e-7, abstol = 1.0e-9, controller = predictive_disco_controller(RadauIIA5()))
#  145.542 μs (1641 allocations: 64.64 KiB)
sol_no_disco_radau = solve(prob_multi, RadauIIA5(); callback = cb_multi2, reltol = 1.0e-7, abstol = 1.0e-9)
#  125.667 μs (1044 allocations: 46.12 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback = cb_multi, reltol = 1.0e-7, abstol = 1.0e-9, controller = PI_disco_controller(Rodas5P()))
#  253.291 μs (1758 allocations: 75.62 KiB)
sol_no_disco_rosenbrock = solve(prob_multi, Rodas5P(); callback = cb_multi2, reltol = 1.0e-7, abstol = 1.0e-9)
#  239.292 μs (947 allocations: 43.89 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_multi, Tsit5(); callback = cb_multi, reltol = 1.0e-7, abstol = 1.0e-9, controller = PI_disco_controller(Tsit5()))
#  99.250 μs (1603 allocations: 59.85 KiB)
sol_no_disco_tsit5 = solve(prob_multi, Tsit5(); callback = cb_multi2, reltol = 1.0e-7, abstol = 1.0e-9)
#  88.000 μs (1105 allocations: 45.27 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

#TEST 4: STIFF DISCONTINUITY
# very stiff discontinuous system
function f_stiff_disc!(du, u, p, t)
    λ = p[1] # stiffness parameter
    return if u[1] < 0.5
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
sol_disco_radau = solve(prob_stiff, RadauIIA5(); callback = cb_stiff, reltol = 1.0e-9, abstol = 1.0e-11, controller = predictive_disco_controller(RadauIIA5()))
#  119.250 μs (1666 allocations: 66.52 KiB)
sol_no_disco_radau = solve(prob_stiff, RadauIIA5(); callback = cb_stiff_f, reltol = 1.0e-9, abstol = 1.0e-11)
#  124.167 μs (1476 allocations: 61.16 KiB)
@test sol_disco_radau.retcode == ReturnCode.Success
@test sol_disco_radau.stats.nreject <= sol_no_disco_radau.stats.nreject

sol_disco_tsit5 = solve(prob_stiff, Tsit5(); callback = cb_stiff, reltol = 1.0e-9, abstol = 1.0e-11, controller = PI_disco_controller(Tsit5()))
#  80.375 μs (1870 allocations: 70.51 KiB)
sol_no_disco_tsit5 = solve(prob_stiff, Tsit5(); callback = cb_stiff_f, reltol = 1.0e-9, abstol = 1.0e-11)
#  79.583 μs (1783 allocations: 68.51 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject

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
cb_daef = ContinuousCallback(cond_dae, default_affect!; maybe_discontinuity = false)

radau_disco = solve(prob_dae, RadauIIA5(); callback = cb_dae, reltol = 1.0e-8, abstol = 1.0e-10, controller = predictive_disco_controller(RadauIIA5()))
#  62.875 μs (754 allocations: 33.98 KiB)
radau_no_disco = solve(prob_dae, RadauIIA5(); callback = cb_daef, reltol = 1.0e-8, abstol = 1.0e-10)
#  58.417 μs (639 allocations: 30.61 KiB)
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
cb_vec2 = VectorContinuousCallback(condition_vec!, default_affect_vec!, 2; maybe_discontinuity = false)

sol_disco_rosenbrock = solve(prob_vec, Rodas5P(); callback = cb_vec, reltol = 1.0e-7, abstol = 1.0e-9, controller = PI_disco_controller(Rodas5P()))
#  175.292 μs (1213 allocations: 54.08 KiB)
sol_no_disco_rosenbrock = solve(prob_vec, Rodas5P(); callback = cb_vec2, reltol = 1.0e-7, abstol = 1.0e-9)
#  194.709 μs (888 allocations: 46.53 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec, reltol = 1.0e-7, abstol = 1.0e-9, controller = PI_disco_controller(Tsit5()))
#  73.333 μs (1349 allocations: 57.71 KiB)
sol_no_disco_tsit5 = solve(prob_vec, Tsit5(); callback = cb_vec2, reltol = 1.0e-7, abstol = 1.0e-9)
#  74.708 μs (1256 allocations: 58.48 KiB)
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
cb2 = ContinuousCallback(cond, default_affect!; maybe_discontinuity = false)

sol_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb, reltol = 1.0e-8, abstol = 1.0e-10, controller = PI_disco_controller(Rodas5P()))
#  215.792 μs (1621 allocations: 61.61 KiB)
sol_no_disco_rosenbrock = solve(prob, Rodas5P(); callback = cb2, reltol = 1.0e-8, abstol = 1.0e-10)
#  172.333 μs (845 allocations: 43.62 KiB)
@test sol_disco_rosenbrock.retcode == ReturnCode.Success
@test sol_disco_rosenbrock.stats.nreject <= sol_no_disco_rosenbrock.stats.nreject

sol_disco_tsit5 = solve(prob, Tsit5(); callback = cb, reltol = 1.0e-8, abstol = 1.0e-10, controller = PI_disco_controller(Tsit5()))
#  66.292 μs (1450 allocations: 59.10 KiB)
sol_no_disco_tsit5 = solve(prob, Tsit5(); callback = cb2, reltol = 1.0e-8, abstol = 1.0e-10)
#  52.500 μs (1117 allocations: 52.55 KiB)
@test sol_disco_tsit5.retcode == ReturnCode.Success
@test sol_disco_tsit5.stats.nreject <= sol_no_disco_tsit5.stats.nreject
