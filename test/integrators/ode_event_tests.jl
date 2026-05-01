using OrdinaryDiffEq,
    RecursiveArrayTools, Test, StaticArrays, DiffEqCallbacks,
    SparseArrays,
    LinearAlgebra
using OrdinaryDiffEqLinear, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqRosenbrock

f = function (u, p, t)
    return -u + sin(-t)
end

prob = ODEProblem(f, 1.0, (0.0, -10.0))

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return -t - 2.95
end

affect! = function (integrator)
    return integrator.u = integrator.u + 2
end

callback = ContinuousCallback(condition, affect!)

sol = solve(prob, Tsit5(), callback = callback)
@test length(sol.t) < 20

# Force integrator to step on event
sol = solve(prob, Tsit5(), callback = callback, tstops = [-2.95])
@test sol(-2.95, continuity = :right) ≈ sol(-2.95, continuity = :left) + 2

condition = function (out, u, t, integrator) # Event when event_f(u,t,k) == 0
    return out[1] = -t - 2.95
end

affect! = function (integrator, events)
    for (idx, dir) in enumerate(events)
        iszero(dir) && continue
        if idx == 1
            integrator.u = integrator.u + 2
        end
    end
    return
end

callback = VectorContinuousCallback(condition, affect!, 1)

sol = solve(prob, Tsit5(), callback = callback)

# Force integrator to step on event
sol = solve(prob, Tsit5(), callback = callback, tstops = [-2.95])
@test sol(-2.95, continuity = :right) ≈ sol(-2.95, continuity = :left) + 2

f = function (du, u, p, t)
    return du[1] = -u[1] + sin(t)
end

prob = ODEProblem(f, [1.0], (0.0, 10.0))

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return t - 2.95
end

affect! = function (integrator)
    return integrator.u = integrator.u .+ 2
end

callback = ContinuousCallback(condition, affect!)

sol = solve(prob, Tsit5(), callback = callback, abstol = 1.0e-8, reltol = 1.0e-6)

# Force integrator to step on event
sol = solve(prob, Tsit5(), callback = callback, abstol = 1.0e-8, reltol = 1.0e-6, tstops = [2.95])
@test sol(2.95, continuity = :right)[1] ≈ sol(2.95, continuity = :left)[1] + 2

condition = function (out, u, t, integrator) # Event when event_f(u,t,k) == 0
    return out[1] = t - 2.95
end

affect! = function (integrator, events)
    for (idx, dir) in enumerate(events)
        iszero(dir) && continue
        if idx == 1
            integrator.u = integrator.u .+ 2
        end
    end
    return
end

callback = VectorContinuousCallback(condition, affect!, 1)

sol = solve(prob, Tsit5(), callback = callback, abstol = 1.0e-8, reltol = 1.0e-6)

# Force integrator to step on event
sol = solve(prob, Tsit5(), callback = callback, abstol = 1.0e-8, reltol = 1.0e-6, tstops = [2.95])
@test sol(2.95, continuity = :right)[1] ≈ sol(2.95, continuity = :left)[1] + 2

f = function (du, u, p, t)
    du[1] = u[2]
    return du[2] = -9.81
end

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return u[1]
end

affect! = nothing
affect_neg! = function (integrator)
    return integrator.u[2] = -integrator.u[2]
end

callback = ContinuousCallback(condition, affect!, affect_neg!, interp_points = 100)

u0 = [50.0, 0.0]
tspan = (0.0, 15.0)
prob = ODEProblem(f, u0, tspan)

sol = solve(prob, Tsit5(), callback = callback, adaptive = false, dt = 1 / 4)

condition = function (out, u, t, integrator) # Event when event_f(u,t,k) == 0
    return out[1] = u[1]
end

affect! = function (integrator, events)
    for (idx, dir) in enumerate(events)
        if dir == 1 && idx == 1 # downcrossing
            integrator.u[2] = -integrator.u[2]
        end
    end
    return
end
affect_neg! = function (integrator)
    return integrator.u[2] = -integrator.u[2]
end

vcb = VectorContinuousCallback(condition, affect!, 1, interp_points = 100)

sol = solve(prob, Tsit5(), callback = vcb, adaptive = false, dt = 1 / 4)

condition_single = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return u
end

affect! = nothing
affect_neg! = function (integrator)
    return integrator.u[2] = -integrator.u[2]
end

callback_single = ContinuousCallback(
    condition_single, affect!, affect_neg!,
    interp_points = 100, idxs = 1
)

u0 = [50.0, 0.0]
tspan = (0.0, 15.0)
prob = ODEProblem(f, u0, tspan)

sol = solve(prob, Tsit5(), callback = callback_single, adaptive = false, dt = 1 / 4)
sol = solve(prob, Tsit5(), callback = callback_single, save_everystep = false)
t = sol.t[end ÷ 2] # this is the callback time point
sol = solve(prob, Tsit5(), callback = callback_single, saveat = t)
@test count(x -> x == t, sol.t) == 2
sol = solve(prob, Tsit5(), callback = callback_single, saveat = t - eps(t))
@test count(x -> x == t, sol.t) == 2
# check interpolation @ discontinuity
@test sol(t, continuity = :right)[2] > 0
@test sol(t, continuity = :left)[2] < 0

#plot(sol,denseplot=true)

condition_single = function (out, u, t, integrator) # Event when event_f(u,t,k) == 0
    return out[1] = u[1]
end

affect! = function (integrator, events)
    for (idx, dir) in enumerate(events)
        if dir == 1 && idx == 1 # downcrossing
            integrator.u[2] = -integrator.u[2]
        end
    end
    return
end

callback_single = VectorContinuousCallback(
    condition_single, affect!, 1,
    interp_points = 100
)

sol = solve(prob, Tsit5(), callback = callback_single, adaptive = false, dt = 1 / 4)
sol = solve(prob, Tsit5(), callback = callback_single, save_everystep = false)
t = sol.t[end ÷ 2] # this is the callback time point
sol = solve(prob, Tsit5(), callback = callback_single, saveat = t)
@test count(x -> x == t, sol.t) == 2
sol = solve(prob, Tsit5(), callback = callback_single, saveat = t - eps(t))
@test count(x -> x == t, sol.t) == 2
# check interpolation @ discontinuity
@test sol(t, continuity = :right)[2] > 0
@test sol(t, continuity = :left)[2] < 0

sol = solve(prob, Vern6(), callback = callback)
sol = solve(prob, Vern6(), callback = vcb)
#plot(sol,denseplot=true)
sol = solve(prob, BS3(), callback = vcb)
sol = solve(prob, BS3(), callback = callback)

sol33 = solve(prob, Vern7(), callback = callback)
sol33 = solve(prob, Vern7(), callback = vcb)

bounced = ODEProblem(f, sol.u[end - 1], (0.0, 1.0))
sol_bounced = solve(bounced, Vern6(), callback = callback, dt = sol.t[end] - sol.t[end - 1])
#plot(sol_bounced,denseplot=true)
sol_bounced(0.04) # Complete density
@test maximum(
    maximum.(
        map(
            (i) -> sol.k[end][i] - sol_bounced.k[2][i],
            1:length(sol.k[end])
        )
    )
) ==
    0

sol2 = solve(prob, Vern6(), callback = callback, adaptive = false, dt = 1 / 2^4)
#plot(sol2)

sol2 = solve(prob, Vern6())

sol3 = solve(prob, Vern6(), saveat = [0.5])

## Saving callback

condition = function (u, t, integrator)
    return true
end
affect! = function (integrator) end

save_positions = (true, false)
saving_callback = DiscreteCallback(condition, affect!, save_positions = save_positions)

sol4 = solve(prob, Tsit5(), callback = saving_callback)

@test sol2(3) ≈ sol(3)

affect! = function (integrator)
    return derivative_discontinuity!(integrator, false)
end
saving_callback2 = DiscreteCallback(condition, affect!, save_positions = save_positions)
sol4 = solve(prob, Tsit5(), callback = saving_callback2)

cbs = CallbackSet(saving_callback, saving_callback2)
sol4_extra = solve(prob, Tsit5(), callback = cbs)

@test length(sol4_extra.t) == 2length(sol4.t) - 1

condition = function (u, t, integrator)
    return u[1]
end

vcondition = function (out, u, t, integrator)
    return out[1] = u[1]
end

affect! = function (integrator, retcode = nothing)
    return if retcode === nothing
        terminate!(integrator)
    else
        terminate!(integrator, retcode)
    end
end

vaffect! = function (integrator, events, retcode = nothing)
    for (idx, dir) in enumerate(events)
        iszero(dir) && continue
        if idx == 1
            if retcode === nothing
                terminate!(integrator)
            else
                terminate!(integrator, retcode)
            end
        end
    end
    return
end

terminate_callback = ContinuousCallback(condition, affect!)
custom_retcode_callback = ContinuousCallback(
    condition,
    x -> affect!(x, ReturnCode.MaxIters)
)
vterminate_callback = VectorContinuousCallback(vcondition, vaffect!, 1)
vcustom_retcode_callback = VectorContinuousCallback(
    vcondition,
    (x, events) -> vaffect!(
        x, events,
        ReturnCode.MaxIters
    ),
    1
)

tspan2 = (0.0, Inf)
prob2 = ODEProblem(f, u0, tspan2)

sol5 = solve(prob2, Tsit5(), callback = terminate_callback)
sol5_1 = solve(prob2, Tsit5(), callback = custom_retcode_callback)

@test sol5.retcode == ReturnCode.Terminated
@test sol5_1.retcode == ReturnCode.MaxIters
@test sol5.u[end][1] < 3.0e-12
@test sol5.t[end] ≈ sqrt(50 * 2 / 9.81)

sol5 = solve(prob2, Tsit5(), callback = vterminate_callback)
sol5_1 = solve(prob2, Tsit5(), callback = vcustom_retcode_callback)

@test sol5.retcode == ReturnCode.Terminated
@test sol5_1.retcode == ReturnCode.MaxIters
@test sol5.u[end][1] < 3.0e-12
@test sol5.t[end] ≈ sqrt(50 * 2 / 9.81)

affect2! = function (integrator)
    return if integrator.t >= 3.5
        terminate!(integrator)
    else
        integrator.u[2] = -integrator.u[2]
    end
end

vaffect2! = function (integrator, events)
    for (idx, dir) in enumerate(events)
        if dir == 1 && idx == 1 # downcrossing
            if integrator.t >= 3.5
                terminate!(integrator)
            else
                integrator.u[2] = -integrator.u[2]
            end
        end
    end
    return
end

terminate_callback2 = ContinuousCallback(condition, nothing, affect2!, interp_points = 100)
vterminate_callback2 = VectorContinuousCallback(
    vcondition, vaffect2!, 1,
    interp_points = 100
)

sol5 = solve(prob2, Vern7(), callback = terminate_callback2)

@test sol5.u[end][1] < 1.3e-10
@test sol5.t[end] ≈ 3 * sqrt(50 * 2 / 9.81)

sol5 = solve(prob2, Vern7(), callback = vterminate_callback2)

@test sol5.u[end][1] < 1.3e-10
@test sol5.t[end] ≈ 3 * sqrt(50 * 2 / 9.81)

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return t - 4
end

affect! = function (integrator)
    return terminate!(integrator)
end

terminate_callback3 = ContinuousCallback(condition, affect!, interp_points = 1000)

bounce_then_exit = CallbackSet(callback, terminate_callback3)

sol6 = solve(prob2, Vern7(), callback = bounce_then_exit)

@test sol6.u[end][1] > 0
@test sol6.u[end][1] < 100
@test sol6.t[end] ≈ 4

# Test ContinuousCallback hits values on the steps
t_event = 100.0
f_simple(u, p, t) = 1.00001 * u
event_triggered = false
condition_simple(u, t, integrator) = t_event - t
function affect_simple!(integrator)
    global event_triggered
    return event_triggered = true
end
cb = ContinuousCallback(condition_simple, nothing, affect_simple!)
prob = ODEProblem(f_simple, [1.0], (0.0, 2.0 * t_event))
sol = solve(prob, Tsit5(), callback = cb, adaptive = false, dt = 10.0)
@test event_triggered

# https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/328
ode = ODEProblem((du, u, p, t) -> (@. du .= -u), ones(5), (0.0, 100.0))
sol = solve(ode, AutoTsit5(Rosenbrock23()), callback = TerminateSteadyState())
sol1 = solve(ode, Tsit5(), callback = TerminateSteadyState())
@test sol.u == sol1.u

# DiscreteCallback
f = function (du, u, p, t)
    du[1] = -0.5 * u[1] + 10
    return du[2] = -0.5 * u[2]
end

u0 = [10, 10.0]
tstop = [5.0; 8.0]
prob = ODEProblem(f, u0, (0, 10.0))
condition = (u, t, integrator) -> t in tstop
affect! = (integrator) -> integrator.u .= 1.0
save_positions = (true, true)
cb = DiscreteCallback(condition, affect!, save_positions = save_positions)
sol1 = solve(prob, Tsit5(), callback = cb, tstops = tstop, saveat = tstop)
@test count(x -> x == tstop[1], sol1.t) == 2
@test count(x -> x == tstop[2], sol1.t) == 2
sol2 = solve(prob, Tsit5(), callback = cb, tstops = tstop, saveat = prevfloat.(tstop))
@test count(x -> x == tstop[1], sol2.t) == 2
@test count(x -> x == tstop[2], sol2.t) == 2

# check VectorContinuousCallback works for complex valued solutions
# see issue #1222
f = (u, p, t) -> -1.0im * u
prob = ODEProblem(f, complex([1.0]), (0.0, 1.0))
condition = function (out, u, t, integrator)
    return out[1] = t - 0.5
end
n = 0
affect! = function (integrator, events)
    global n
    for (_, dir) in enumerate(events)
        iszero(dir) && continue
        n += 1
    end
    return
end
callback = VectorContinuousCallback(condition, affect!, 1)
sol = solve(prob, Tsit5(), callback = callback)
@test n == 1

# case of immutable partitioned state
f = function (u, p, t)
    return ArrayPartition(SVector{1}(u[2]), SVector{1}(-9.81))
end

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    return u[1]
end

affect! = nothingf = affect_neg! = function (integrator)
    return integrator.u = ArrayPartition(SVector{1}(integrator.u[1]), SVector{1}(-integrator.u[2]))
end

callback = ContinuousCallback(condition, affect!, affect_neg!, interp_points = 100)

u0 = ArrayPartition(SVector{1}(50.0), SVector{1}(0.0))
tspan = (0.0, 15.0)
prob = ODEProblem(f, u0, tspan)

sol = solve(prob, Tsit5(), callback = callback, adaptive = false, dt = 1 / 4)

# check that multiple discrete callbacks with save_everystep do not double save
# https://github.com/SciML/DifferentialEquations.jl/issues/711

prob = ODEProblem((u, p, t) -> 1.01u, [1.0], (0.0, 10.0))
savecond1(u, t, integrator) = t == 2.5
savecond2(u, t, integrator) = t == 2.5
saveaffect1!(integrator) = nothing
saveaffect2!(integrator) = nothing
cb1 = DiscreteCallback(savecond1, saveaffect1!, save_positions = (false, false))
cb2 = DiscreteCallback(savecond2, saveaffect2!, save_positions = (false, false))
cb = CallbackSet(cb1, cb2)
sol = solve(prob, Tsit5(), callback = cb, tstops = [2.5])
@test !any(diff(sol.t) .== 0)

# DifferentialEquations 848
prob = ODEProblem((x, p, t) -> -1.01 * x, ones(2), (0.0, 1.0))
integrator = init(prob, Tsit5(), save_everystep = false)
set_u!(integrator, 2 * ones(2))
step!(integrator, 1.0e-5, true)
@test all(u -> u > 1.5, integrator.u)

# https://github.com/SciML/OrdinaryDiffEq.jl/pull/1777
@testset "Callbacks with LinearExponential" begin
    A = Matrix(sprand(ComplexF64, 100, 100, 0.5))
    A += A'

    t_l = LinRange(0, 1, 100)

    saved_values = SavedValues(Float64, Float64)
    function save_func(u, t, integrator)
        real(u' * A * u)
    end
    cb = SavingCallback(save_func, saved_values, saveat = t_l)

    u0 = normalize(rand(ComplexF64, 100))
    A = MatrixOperator(-1im * A)
    prob = ODEProblem(A, u0, (0, 1.0))
    solve(prob, LinearExponential(), dt = t_l[2] - t_l[1], callback = cb)
    @test length(saved_values.saveval) == length(t_l)
end

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/3222
# Friction state machine where a velocity zero crossing under one discrete state
# (SLIDING_LEFT, condition 6) transitions to a new state (SLIDING_RIGHT) where
# the corresponding condition (condition 2) is immediately at zero. The rootfinder
# can't detect this because condition 2 starts at zero rather than crossing through
# it. With the simultaneous events API, the user can chain the transition in the
# affect! function.
@testset "Issue #3222: chained velocity crossing in friction state machine" begin
    _LEFT_STOP, _SLIDING_RIGHT, _STOPPED, _SLIDING_LEFT, _RIGHT_STOP = 1:5

    mutable struct _FrictionParams
        M::Float64
        kf::Float64
        Fsp::Float64
        xls::Float64
        xrs::Float64
        Fsf::Float64
        Fkf::Float64
        cor::Float64
        state::Int
        A::Float64
        minbouncespeed::Float64
    end

    _P(t) = t > 0.0 ? 6895 * (18 + 4 * sin(2π * 250 * t)) : 0.0

    function _friction_ode(u, p, t)
        v, x = u
        Fspring = p.kf * x + p.Fsp
        st = p.state
        if st == _LEFT_STOP
            dv = 0.0; dx = 0.0
        elseif st == _SLIDING_RIGHT
            dv = 1 / p.M * (_P(t) * p.A - p.Fkf - Fspring); dx = v
        elseif st == _STOPPED
            dv = 0.0; dx = 0.0
        elseif st == _SLIDING_LEFT
            dv = 1 / p.M * (_P(t) * p.A + p.Fkf - Fspring); dx = v
        else
            dv = 0.0; dx = 0.0
        end
        @SVector [dv, dx]
    end

    function _friction_condition(out, u, t, integ)
        p = integ.p
        v, x = u
        Fspring = p.kf * x + p.Fsp
        st = p.state
        out[1] = (st == _LEFT_STOP) * (_P(t) * p.A - p.Fsf - Fspring)
        out[2] = (st == _SLIDING_RIGHT) * (-v)
        out[3] = (st == _SLIDING_RIGHT) * (x - p.xrs)
        out[4] = (st == _STOPPED) * (_P(t) * p.A - p.Fsf - Fspring)
        out[5] = (st == _STOPPED) * (-(_P(t) * p.A + p.Fsf - Fspring))
        out[6] = (st == _SLIDING_LEFT) * v
        out[7] = (st == _SLIDING_LEFT) * (-(x - p.xls))
        out[8] = (st == _RIGHT_STOP) * (-(_P(t) * p.A - Fspring + p.Fsf))
    end

    affect_log = Tuple{Float64, Int, Int}[]  # (t, idx, new_state)

    function _friction_affect!(integ, events)
        p = integ.p
        v, x = integ.u
        t = integ.t
        Fspring = p.Fsp + x * p.kf

        for (idx, dir) in enumerate(events)
            iszero(dir) && continue

            if idx == 1
                p.state = _SLIDING_RIGHT
            elseif idx == 2
                if _P(t) * p.A - Fspring < -p.Fkf
                    p.state = _SLIDING_LEFT
                else
                    p.state = _STOPPED
                end
            elseif idx == 3
                if p.cor == 0.0 || p.cor * v < p.minbouncespeed
                    p.state = _RIGHT_STOP
                    integ.u = @SVector [0.0, integ.u[2]]
                else
                    p.state = _SLIDING_LEFT
                    integ.u = @SVector [-p.cor * v, integ.u[2]]
                end
            elseif idx == 4
                p.state = _SLIDING_RIGHT
            elseif idx == 5
                p.state = _SLIDING_LEFT
            elseif idx == 6
                if _P(t) * p.A - Fspring > p.Fkf
                    p.state = _SLIDING_RIGHT
                else
                    p.state = _STOPPED
                    integ.u = @SVector [0.0, integ.u[2]]
                end
            elseif idx == 7
                if p.cor == 0.0 || -p.cor * v < p.minbouncespeed
                    p.state = _LEFT_STOP
                    integ.u = @SVector [0.0, integ.u[2]]
                else
                    p.state = _SLIDING_RIGHT
                    integ.u = @SVector [-p.cor * v, integ.u[2]]
                end
            elseif idx == 8
                p.state = _SLIDING_LEFT
            end
            push!(affect_log, (t, idx, p.state))

            # Chain: when a velocity-zero event (idx 2 or 6) transitions to
            # the opposite sliding state, v is still ≈0 so the new state's
            # velocity-crossing condition is immediately at its trigger point.
            # The rootfinder can't detect a crossing that starts at zero, so
            # handle it explicitly here.
            v_now = integ.u[1]
            if (idx == 2 || idx == 6) && abs(v_now) < 1.0e-8 &&
                    (p.state == _SLIDING_RIGHT || p.state == _SLIDING_LEFT)
                prev_state = p.state
                if p.state == _SLIDING_RIGHT
                    if _P(t) * p.A - Fspring < -p.Fkf
                        p.state = _SLIDING_LEFT
                    else
                        p.state = _STOPPED
                    end
                elseif p.state == _SLIDING_LEFT
                    if _P(t) * p.A - Fspring > p.Fkf
                        p.state = _SLIDING_RIGHT
                    else
                        p.state = _STOPPED
                        integ.u = @SVector [0.0, integ.u[2]]
                    end
                end
                push!(affect_log, (t, idx, p.state))
            end
        end
    end

    vcb = VectorContinuousCallback(
        _friction_condition, _friction_affect!, nothing,
        save_positions = (true, true), 8
    )

    p = _FrictionParams(
        0.0003361185937456267, 1488.578099595049, 8.54503372291542,
        0.0, 0.001016, 1.0, 1.0, 0.3, _LEFT_STOP,
        7.373657906574696e-5, 0.001
    )
    u0 = @SVector [0.0, 0.0]
    prob = ODEProblem(_friction_ode, u0, (-10.0e-6, 0.0025), p)
    empty!(affect_log)
    sol = solve(
        prob, Rodas5(), callback = vcb, reltol = 1.0e-6, abstol = 1.0e-6,
        save_everystep = true
    )

    # The chained logic should catch the v≈0 condition when idx=6 transitions
    # to SLIDING_RIGHT, immediately chaining to STOPPED.
    has_chained_transition = any(affect_log) do (t, idx, new_state)
        idx == 6 && new_state == _STOPPED
    end
    @test has_chained_transition

    # With the chain, velocity should stay ≈0 through the previously-problematic
    # region (t=0.00128 to t=0.00135) instead of going negative without a callback.
    @test abs(sol(0.00128)[1]) < 0.01
    @test abs(sol(0.00135)[1]) < 0.01
end
