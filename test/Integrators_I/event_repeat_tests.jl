using OrdinaryDiffEq, Test

f(u, p, t) = u
prob = ODEProblem(f, 1.0, (0.0, 100eps(1.0)))
dt = 4.7eps(1.0)
cond(u, p, t) = u - exp(70eps(1.0))

c = Ref(0)
function affect!(integrator)
    c[] += 1
    return
end
cb = ContinuousCallback(cond, affect!)
@info "Event Repeat Test 1"
sol = solve(prob, Tsit5(); adaptive = false, callback = cb, dt, save_everystep = false)
@test c[] == 1

function condition_v(out, u, t, integrator)
    out[1] = u - exp(50eps(1.0))
    out[2] = u - exp(60eps(1.0))
    out[3] = u - exp(70eps(1.0))
    out[4] = u - exp(80eps(1.0))
    out[5] = u - exp(90eps(1.0))
    return
end

c1 = Ref(0)
c2 = Ref(0)
c3 = Ref(0)
c4 = Ref(0)
c5 = Ref(0)

function affect_v!(integrator, events)
    for (idx, dir) in enumerate(events)
        iszero(dir) && continue
        if idx == 1
            c1[] += 1
        elseif idx == 2
            c2[] += 1
        elseif idx == 3
            c3[] += 1
        elseif idx == 4
            c4[] += 1
        elseif idx == 5
            c5[] += 1
        end
    end
    return
end

cb = VectorContinuousCallback(condition_v, affect_v!, 5)
@info "Event Repeat Test 2"
sol = solve(prob, Tsit5(); adaptive = false, callback = cb, dt, save_everystep = false)

@test c1[] == 1
@test c2[] == 1
@test c3[] == 1
@test c4[] == 1
@test c5[] == 1

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/4991
# A bounce whose flight is shorter than `dt * repeat_nudge` must not be skipped.
function bounce!(du, u, p, t)
    du[1] = u[2]
    du[2] = -9.81
    du[3:end] .= 0
    return
end
@testset "Short flights near the Zeno time, $(nameof(typeof(alg))), e = $e" for alg in (Tsit5(), Vern9()),
        (e, h0) in ((0.8, 1.0), (0.7, 10.0))

    tz = sqrt(2h0 / 9.81) * (1 + e) / (1 - e)
    bounces = Ref(0)
    bounce_cb = ContinuousCallback(
        (u, t, integrator) -> u[1], nothing,
        integrator -> (bounces[] += 1; integrator.u[2] *= -e)
    )
    sol = solve(
        ODEProblem(bounce!, [h0, 0.0], (0.0, tz + 1)), alg;
        callback = bounce_cb, abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test bounces[] > 40
    @test minimum(first, sol.(range(0, tz - 1.0e-3, length = 10^4))) > -1.0e-8
end

# A large state component that the condition does not read must not hide short flights
@testset "Short flights near the Zeno time with a large unrelated state" begin
    tz = sqrt(2 / 9.81) * 9
    bounce_cb = ContinuousCallback(
        (u, t, integrator) -> u[1], nothing, integrator -> (integrator.u[2] *= -0.8)
    )
    sol = solve(
        ODEProblem(bounce!, [1.0, 0.0, 1.0e10], (0.0, tz + 1)), Tsit5();
        callback = bounce_cb, abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test minimum(first, sol.(range(0, tz - 1.0e-3, length = 10^4))) > -1.0e-8
end

# Rounding noise of a large state must not read as new crossings after the Zeno time
@testset "No chatter after the Zeno time with a floor at 1000" begin
    tz = sqrt(2 / 9.81) * 9
    times = Float64[]
    floor_cb = ContinuousCallback(
        (u, t, integrator) -> u[1] - 1000, nothing,
        integrator -> (push!(times, integrator.t); integrator.u[2] *= -0.8)
    )
    solve(
        ODEProblem(bounce!, [1001.0, 0.0], (0.0, tz + 1)), Tsit5();
        callback = floor_cb, abstol = 1.0e-10, reltol = 1.0e-10, save_everystep = false
    )
    @test count(<(tz), times) > 40
    @test count(>=(tz), times) < 100
end

# Components that fired together are each checked, whatever their order
@testset "Simultaneous VectorContinuousCallback components, order $order" for order in ([1, 2], [2, 1])
    function balls!(du, u, p, t)
        du[1] = u[2]
        du[2] = -9.81
        du[3] = u[4]
        du[4] = -9.81
        return
    end
    e = (0.8, 0.005)
    t1 = sqrt(2 / 9.81)
    impacts = [[t1, t1 * (1 + 2e[ball])] for ball in 1:2]
    prob = ODEProblem(balls!, [1.0, 0.0, 1.0, 0.0], (0.0, 1.2))
    vcc_times = [Float64[], Float64[]]
    function vcc_affect!(integrator, events)
        for idx in 1:2
            ball = order[idx]
            if events[idx] != 0 && integrator.u[2ball] < 0
                push!(vcc_times[ball], integrator.t)
                integrator.u[2ball] *= -e[ball]
            end
        end
        return
    end
    vcc = VectorContinuousCallback(
        (out, u, t, integrator) -> (out .= (u[1], u[3])[order]), vcc_affect!, 2
    )
    solve(prob, Tsit5(); callback = vcc, abstol = 1.0e-10, reltol = 1.0e-10)
    set_times = [Float64[], Float64[]]
    ball_cb(ball) = ContinuousCallback(
        (u, t, integrator) -> u[2ball - 1], nothing,
        integrator -> (push!(set_times[ball], integrator.t); integrator.u[2ball] *= -e[ball])
    )
    solve(prob, Tsit5(); callback = CallbackSet(ball_cb(1), ball_cb(2)), abstol = 1.0e-10, reltol = 1.0e-10)
    for ball in 1:2
        @test vcc_times[ball][1:2] ≈ impacts[ball] rtol = 1.0e-6
        @test set_times[ball][1:2] ≈ impacts[ball] rtol = 1.0e-6
    end
end

function still!(du, u, p, t)
    du .= 0
    return
end

# A later crossing inside the nudge window must not be lost
@testset "Three roots in one nudge window, vector = $vector, interp_points = $ip" for vector in (false, true),
        ip in (0, 10)

    times = Float64[]
    g(t) = (t - 0.5) * (t - 0.502) * (t - 0.506)
    cb = if vector
        VectorContinuousCallback(
            (out, u, t, integrator) -> (out[1] = g(t)),
            (integrator, idx) -> push!(times, integrator.t), 1; interp_points = ip
        )
    else
        ContinuousCallback(
            (u, t, integrator) -> g(t), integrator -> push!(times, integrator.t);
            interp_points = ip
        )
    end
    solve(ODEProblem(still!, [0.0], (0.0, 1.0)), Tsit5(); dt = 0.5, adaptive = false, callback = cb)
    @test times ≈ [0.5, 0.502, 0.506]
end

# The nudge search must stay bounded after an exact root at `t = 0`
@testset "Exact BigFloat root at zero, vector = $vector" for vector in (false, true)
    count = Ref(0)
    cb = if vector
        VectorContinuousCallback(
            (out, u, t, integrator) -> (out[1] = t), (integrator, idx) -> (count[] += 1), 1;
            abstol = 0
        )
    else
        ContinuousCallback((u, t, integrator) -> t, integrator -> (count[] += 1); abstol = 0)
    end
    sol = solve(
        ODEProblem(still!, BigFloat[0], (BigFloat(-0.5), BigFloat(0.5))), Tsit5();
        dt = BigFloat(0.5), adaptive = false, callback = cb
    )
    @test sol.retcode == ReturnCode.Success
    @test count[] == 1
    @test sol.stats.ncondition < 200
end

# A short flight seen above the callback tolerance at a single sample must not be skipped
@testset "Short flight with one resolved sample, vector = $vector" for vector in (false, true)
    function ball!(du, u, p, t)
        du[1] = u[2]
        du[2] = -2
        return
    end
    times = Float64[]
    affect!(integrator) = (push!(times, integrator.t); integrator.u[2] *= -0.00045)
    cb = if vector
        VectorContinuousCallback(
            (out, u, t, integrator) -> (out[1] = u[1]),
            (integrator, idx) -> integrator.u[2] < 0 && affect!(integrator), 1; abstol = 5.0e-9
        )
    else
        ContinuousCallback((u, t, integrator) -> u[1], nothing, affect!; abstol = 5.0e-9)
    end
    solve(ODEProblem(ball!, [0.25, 0.0], (0.0, 0.6)), Tsit5(); dt = 0.5, adaptive = false, callback = cb)
    @test times[1:2] ≈ [0.5, 0.50045]
end
