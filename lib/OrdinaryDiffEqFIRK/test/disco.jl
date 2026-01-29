using OrdinaryDiffEqFIRK, DiffEqDevTools, Test, LinearAlgebra

#test example discontinuous at u = 1
f(u, p, t) = u[1] < 1 ? [2u[1]] : [-3u[1] + 5]
u0 = [0.1]           
tspan = (0.0, 1.5)
prob = ODEProblem(f, u0, tspan)

#define callback
condition(u, t, integrator) = u[1] - 1
function affect!(integrator)
    integrator.u[1] += 10
    println("Callback fired at t = ", integrator.t)
end
cb = ContinuousCallback(condition, affect!; is_discontinuity = true)

#disco solve
sol = solve(prob, RadauIIA5(is_disco = true); callback = cb, reltol = 1e-6)
#fixed order solve
sol2 = solve(prob, RadauIIA5(is_disco = false); callback = cb, reltol = 1e-6)

#two discontinuity functions
function f(u, p, t)
    if u[1] < 1
        [2u[1]]                 # region 1: grows to hit u = 1
    elseif u[1] < 2
        [u[1] + 0.2]            # region 2: continues increasing to hit u = 2
    else
        [-4u[1] + 12]           # region 3: after 2, moves toward u ≈ 3
    end
end

u0 = [0.1]
tspan = (0.0, 2.5)
prob = ODEProblem(f, u0, tspan)

#define callbacks
condition1(u, t, integrator) = u[1] - 1
function affect1!(integrator)
    println("Callback 1 fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb1 = ContinuousCallback(condition1, affect1!; is_discontinuity = true)

condition2(u, t, integrator) = u[1] - 2
function affect2!(integrator)
    println("Callback 2 fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb2 = ContinuousCallback(condition2, affect2!; is_discontinuity = true)
cb = CallbackSet(cb1, cb2)

#disco solve
sol = solve(prob, RadauIIA5(is_disco = true); callback = cb, reltol = 1e-6)
#fixed order solve
sol2 = solve(prob, RadauIIA5(is_disco = false); callback = cb, reltol = 1e-6)

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
    println("Multi-exponential discontinuity 1 callback fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb_multi_1 = ContinuousCallback(cond_multi_1, affect_multi_1!; is_discontinuity = true)

cond_multi_2(u, t, integrator) = u[1] - 0.8
function affect_multi_2!(integrator)
    println("Multi-exponential discontinuity 2 callback fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb_multi_2 = ContinuousCallback(cond_multi_2, affect_multi_2!; is_discontinuity = true)
cb_multi = CallbackSet(cb_multi_1, cb_multi_2)

#disco solve
sol_multi_cb_1 = solve(prob_multi, RadauIIA5(is_disco = true); callback=cb_multi, reltol=1e-7, abstol=1e-9)
#fixed order solve
sol = solve(prob_multi, RadauIIA5(is_disco = false); callback=cb_multi, reltol = 1e-7, abstol = 1e-9)

# 2D system with exponential coupling and discontinuity
function f_2d_exp!(du, u, p, t)
    if u[1] + u[2] < 1.0
        du[1] = 2 * exp(u[1]) - u[2]
        du[2] = -3 * u[1] + 4 * exp(u[2])
    else
        du[1] = u[1]
        du[2] = u[2]
    end
end

u0_2d = [0.1, 0.2]
tspan_2d = (0.0, 2.0)
prob_2d = ODEProblem(f_2d_exp!, u0_2d, tspan_2d)

#define callback
cond_2d(u, t, integrator) = u[1] + u[2] - 1.0
function affect_2d!(integrator)
    println("2D exponential discontinuity callback fired at t=$(integrator.t), u=$(integrator.u)")
    @test 0.98 < integrator.u[1] + integrator.u[2] < 1.02
end
cb_2d = ContinuousCallback(cond_2d, affect_2d!; is_discontinuity = true)

#disco solve
sol_2d_cb = solve(prob_2d, RadauIIA5(is_disco = true); callback=cb_2d, reltol=1e-8, abstol=1e-10)
#fixed order solve
sol = solve(prob_2d, RadauIIA5(is_disco = false); callback=cb_2d, reltol = 1e-8, abstol = 1e-10)

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
    println("Stiff discontinuity callback fired at t=$(integrator.t), u=$(integrator.u[1])")
end
cb_stiff = ContinuousCallback(cond_stiff, affect_stiff!; is_discontinuity = true)

#disco solve
sol_stiff_cb = solve(prob_stiff, RadauIIA5(is_disco = true); callback=cb_stiff, reltol=1e-9, abstol=1e-11)
#fixed order solve
sol = solve(prob_stiff, RadauIIA5(is_disco = false); callback=cb_stiff, reltol = 1e-9, abstol = 1e-11)

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
for (i, disc_val) in enumerate(disc_values)
    local cond_func(u, t, integrator) = u[1] - disc_val
    function affect_func!(integrator)
        println("Dense discontinuity $i fired at t=$(integrator.t), u=$(integrator.u[1])")
    end
    push!(cbs_many, ContinuousCallback(cond_func, affect_func!; is_discontinuity = true))
end
cb_many = CallbackSet(cbs_many...)

#disco solve
sol_many_cb = solve(prob_many, RadauIIA5(is_disco = true); callback=cb_many, reltol=1e-10, abstol=1e-12)
#fixed order solve
sol = solve(prob_many, RadauIIA5(is_disco = false); callback=cb_many, reltol=1e-10, abstol=1e-12)

# discontinuity in u (state jump)
function f_state_jump!(du, u, p, t)
    du[1] = 2 * u[1]  # smooth exponential growth until discontinuity
end

u0_jump = [0.5]
tspan_jump = (0.0, 2.0)
prob_jump = ODEProblem(f_state_jump!, u0_jump, tspan_jump)

# define callback that causes a sudden jump in state
cond_jump(u, t, integrator) = u[1] - 1.0
function affect_jump!(integrator)
    integrator.u[1] = 0.2  # sudden jump from ~1.0 down to 0.2
    println("State discontinuity callback fired at t=$(integrator.t), u jumps to 0.2")
end
cb_jump = ContinuousCallback(cond_jump, affect_jump!; is_discontinuity = true)

#disco solve
sol_jump = solve(prob_jump, RadauIIA5(is_disco = true); callback=cb_jump, reltol=1e-8, abstol=1e-10)
#fixed order solve
sol_jump2 = solve(prob_jump, RadauIIA5(is_disco = false); callback=cb_jump, reltol=1e-8, abstol=1e-10)
