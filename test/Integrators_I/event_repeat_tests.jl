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
