using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF, OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner

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

sol_disco = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  277.833 μs (8033 allocations: 251.14 KiB)
sol_no_disco = solve(prob, RadauIIA5(); callback = cb2, reltol = 1e-6)
#  343.041 μs (10008 allocations: 311.02 KiB)

@profview for i in 1:1000 
    solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
end

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
sol_disco = solve(prob, RadauIIA5(); callback = cb, reltol = 1e-6)
#  1.664 ms (43703 allocations: 1.35 MiB)
#fixed order solve
sol_no_disco = solve(prob, RadauIIA5(); callback = cb2, reltol = 1e-6)
#   1.266 ms (37019 allocations: 1.12 MiB)



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
#    202.834 μs (2770 allocations: 93.23 KiB)
#fixed order solve
sol_no_disco = solve(prob_multi, RadauIIA5(); callback=cb_multi2, reltol = 1e-7, abstol = 1e-9)
#  122.875 μs (1136 allocations: 54.52 KiB)

@profview for i in 1:1000 
    solve(prob_multi, RadauIIA5(); callback = cb_multi, reltol = 1e-6)
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
#  131.875 μs (1956 allocations: 74.03 KiB)
#fixed order solve
sol_no_disco = solve(prob_stiff, RadauIIA5(); callback=cb_stiff_f, reltol = 1e-9, abstol = 1e-11)
#  119.417 μs (1480 allocations: 59.55 KiB)


#TEST 5: MULTIPLE DISCONTINUITIES IN SMALL RANGE
# multiple discontinuities in very small range (1e-6 apart, 5 discontinuities)
function f_many_disc!(du, u, p, t)
    du[1] = u[1] + 1  # simple linear growth
end

u0_many = [0.0]
tspan_many = (0.0, 1.0)
prob_many = ODEProblem(f_many_disc!, u0_many, tspan_many)

# create 5 discontinuities spaced 1e-6 apart
disc_values = [0.1 + i * 1e-6 for i = 0:4]

# define callbacks for each discontinuity
cbs_many = []
cbs_many_f = []
for (i, disc_val) in enumerate(disc_values)
    local cond_func(u, t, integrator) = u[1] - disc_val
    function affect_func!(integrator)
        #println("Dense discontinuity $i fired at t=$(integrator.t), u=$(integrator.u[1])")
    end
    push!(cbs_many, ContinuousCallback(cond_func, affect_func!; is_discontinuity = true))
    push!(cbs_many_f, ContinuousCallback(cond_func, affect_func!; is_discontinuity = false))
end
cb_many = CallbackSet(cbs_many...)
cb_many_f = CallbackSet(cbs_many_f...)

#disco solve
sol_disco = solve(prob_many, RadauIIA5(); callback=cb_many, reltol=1e-10, abstol=1e-12)
#   111.333 μs (907 allocations: 36.94 KiB)
#fixed order solve
sol_no_disco = solve(prob_many, RadauIIA5(); callback=cb_many_f, reltol=1e-10, abstol=1e-12)
#   111.666 μs (907 allocations: 36.94 KiB)

#TEST 6: DISCONTINUOUS DAE
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

radau_no_disco = solve(prob_dae, RadauIIA5(); callback=cb_daef, reltol=1e-8, abstol=1e-10)
#  83.500 μs (769 allocations: 35.72 KiB)
radau_disco = solve(prob_dae, RadauIIA5(); callback=cb_dae, reltol=1e-8, abstol=1e-10)
#  101.542 μs (1230 allocations: 48.16 KiB)
 
#TEST 7: VECTOR CALLBACK
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

cb = VectorContinuousCallback(condition!, affect!, 2;) 
cb2 = VectorContinuousCallback(condition!, affect!, 2; is_discontinuity = false) 

sol_disco = solve(prob, RadauIIA5(); callback = cb)
#   62.041 μs (849 allocations: 41.64 KiB)
sol_no_disco = solve(prob, RadauIIA5(); callback = cb2)
#   37.375 μs (531 allocations: 25.23 KiB)