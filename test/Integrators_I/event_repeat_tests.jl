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
    prob = ODEProblem(balls!, [1.0, 0.0, 1.0, 0.0], (0.0, 1.2))
    vcc_bounces = [0, 0]
    function vcc_affect!(integrator, events)
        for idx in 1:2
            ball = order[idx]
            if events[idx] != 0 && integrator.u[2ball] < 0
                vcc_bounces[ball] += 1
                integrator.u[2ball] *= -e[ball]
            end
        end
        return
    end
    vcc = VectorContinuousCallback(
        (out, u, t, integrator) -> (out .= (u[1], u[3])[order]), vcc_affect!, 2
    )
    solve(prob, Tsit5(); callback = vcc, abstol = 1.0e-10, reltol = 1.0e-10)
    set_bounces = [0, 0]
    ball_cb(ball) = ContinuousCallback(
        (u, t, integrator) -> u[2ball - 1], nothing,
        integrator -> (set_bounces[ball] += 1; integrator.u[2ball] *= -e[ball])
    )
    solve(prob, Tsit5(); callback = CallbackSet(ball_cb(1), ball_cb(2)), abstol = 1.0e-10, reltol = 1.0e-10)
    @test vcc_bounces == set_bounces
end
