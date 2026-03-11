using StochasticDiffEq, Test
function f(du, u, p, t)
    du[1] = u[2]
    return du[2] = -9.81
end

function g(du, u, p, t)
    return nothing
end

function condition(u, t, integrator) # Event when event_f(u,p,t,k) == 0
    return u[1]
end

affect! = nothing
function affect_neg!(integrator)
    return integrator.u[2] = -integrator.u[2]
end

# Continuous callback
callback = ContinuousCallback(condition, affect!, affect_neg!)

u0 = [50.0, 0.0]
tspan = (0.0, 15.0)
prob = SDEProblem(f, g, u0, tspan)

sol = solve(prob, SRIW1(), callback = callback, adaptive = false, dt = 3 / 4)

@test minimum([u[1] for u in sol.u]) > -1.0e-12 && minimum([u[1] for u in sol.u]) < 1.0e-12

sol = solve(prob, SRIW1(), callback = callback, save_everystep = false)
t = sol.t[end ÷ 2] # this is the callback time point
sol = solve(prob, SRIW1(), callback = callback, saveat = t)
@test count(x -> x == t, sol.t) == 2
sol = solve(prob, SRIW1(), callback = callback, saveat = t - eps(t))
@test count(x -> x == t, sol.t) == 2

function g(du, u, p, t)
    return du[2] = 0.125 * u[2]
end

prob = SDEProblem(f, g, u0, tspan)

sol = solve(prob, SRIW1(), callback = callback)

sol = solve(prob, EM(), callback = callback, dt = 1 / 4)

# Discrete callback
tstop = [5.0; 8.0]
condition_dc = (u, t, integrator) -> t in tstop
affect!_dc = (integrator) -> integrator.u .= 1.0
save_positions = (true, true)
times_finalize_called = 0
callback_dc = DiscreteCallback(
    condition_dc, affect!_dc, save_positions = save_positions,
    finalize = (args...) -> global times_finalize_called += 1
)
sol = solve(prob, SRIW1(), callback = callback_dc, tstops = tstop, saveat = tstop)
@test count(x -> x == tstop[1], sol.t) == 2
@test count(x -> x == tstop[2], sol.t) == 2
@test times_finalize_called == 1
sol = solve(
    prob, SRIW1(), callback = callback_dc, tstops = tstop, saveat = prevfloat.(tstop)
)
@test count(x -> x == tstop[1], sol.t) == 2
@test count(x -> x == tstop[2], sol.t) == 2
@test times_finalize_called == 2

###
# https://github.com/SciML/DifferentialEquations.jl/issues/802
###

function HM_neuron!(du, u, Params, t)
    # Membrane voltage
    du[1] = -u[1]^3 + 3.0 * u[1]^2 + u[2] - u[3] + u[4] * u[1]

    # Fast channel
    du[2] = 1.0 - 5.0 * u[1]^2 - u[2]

    # Slow channel
    val = Params.s * (u[1] + 8.0 / 5.0)
    du[3] = Params.r * (val - u[3])

    # Synapse
    return du[4] = -u[4] / Params.syntau
end

function HM_noise!(du, u, Params, t)
    du[1] = 0.1
    du[2] = 0.0
    du[3] = 0.0
    return du[4] = 0.0
end

tvals = range(0.0, stop = 1999.9, length = 20000);
struct P
    r::Float64
    s::Float64
    I::Float64
    syntau::Float64
    gSyn::Float64
end

params = P(0.001, 1, 0.0, 5, 0.0);
x0 = Float64[-1.6, -11.8, 0.0, 0.0];

condition2(u, t, integrator) = u[1];
affect2!(integrator) = integrator.u[4] = integrator.u[4] - params.gSyn; # Inhibitory exponential decay synapse
cb = ContinuousCallback(condition2, affect2!, nothing);
prob = SDEProblem(HM_neuron!, HM_noise!, x0, (0.0, 1999.9), params);
sol = solve(prob, ImplicitEM(), reltol = 1.0e-4, abstol = 1.0e-6, dense = true, callback = cb);
sol = solve(prob, SKenCarp(), reltol = 1.0e-4, abstol = 1.0e-6, dense = true, callback = cb);

using DiffEqCallbacks

function f(du, u, p, t)
    return du[1] = p[1] - u[1]
end

function g(du, u, p, t)
    return du[1] = p[2]
end

sprob = SDEProblem(f, g, [1.0], (0.0, 10.0), [1.0, 0.1])
sol = solve(sprob, ImplicitEM(), callback = PositiveDomain())

###
# https://github.com/SciML/DifferentialEquations.jl/issues/1121
# saveat + ContinuousCallback interpolation bug fix
###

using Random

@testset "saveat + ContinuousCallback interpolation (issue #1121)" begin
    α = 1.0
    β = 1.0
    σφ = 0.5

    u0_test = [0.5, 0.0]
    tspan_test = (0.0, 10.0)

    function f_test(u, p, t)
        A, φ = u
        return [α * A * cos(φ), 0.0]
    end

    function g_test(u, p, t)
        A, φ = u
        return [β * A, σφ]
    end

    prob_test = SDEProblem(f_test, g_test, u0_test, tspan_test)

    rng = MersenneTwister(123)
    threshold = Ref(rand(rng, 0.2:0.001:1.0))
    cond_test(u, t, integrator) = u[1] - threshold[]
    function affect_test!(integrator)
        integrator.u[1] *= 2.0
        threshold[] = rand(rng, 0.2:0.001:1.0)
    end

    cb_test = ContinuousCallback(cond_test, affect_test!)

    sol_saveat = solve(
        prob_test, SRIW1(), seed = 123, abstol = 1.0e-2, reltol = 1.0e-2,
        maxiters = Int(1.0e10), saveat = 0.01, callback = cb_test
    )

    # Reset RNG for same callback behavior
    rng = MersenneTwister(123)
    threshold = Ref(rand(rng, 0.2:0.001:1.0))

    sol_dense = solve(
        prob_test, SRIW1(), seed = 123, abstol = 1.0e-2, reltol = 1.0e-2,
        maxiters = Int(1.0e10), callback = cb_test
    )

    # Compare saveat values with dense solution sampled at same times
    phase_saveat = [sol_saveat.u[i][2] for i in 1:length(sol_saveat)]
    phase_dense_sampled = [sol_dense(t)[2] for t in sol_saveat.t]

    # The values should match to machine precision
    max_diff = maximum(abs.(phase_saveat .- phase_dense_sampled))
    @test max_diff < 1.0e-10

    # Also verify no large jumps in phase (the original bug symptom)
    dt = 0.01
    expected_std = σφ * sqrt(dt)
    diffs_saveat = diff(phase_saveat)
    @test all(abs.(diffs_saveat) .< 10 * expected_std)
end

###
# Test that save_discretes defaults to true for SDE integrators,
# matching ODE behavior.
###

@testset "save_discretes defaults to true" begin
    f_sd(u, p, t) = 0.5u
    g_sd(u, p, t) = 0.1u
    sprob = SDEProblem(f_sd, g_sd, 1.0, (0.0, 1.0))
    si = init(sprob, ImplicitEM(); dt = 0.01)
    @test si.opts.save_discretes == true
end
