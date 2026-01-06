using OrdinaryDiffEq, Test

f(u, p, t) = u
prob = ODEProblem(f, 1.0, (0.0, 100eps(1.0)))
dt = 4.7eps(1.0)
cond(u, p, t) = u - exp(70eps(1.0))

c = Ref(0)
function affect!(integrator)
    return c[] += 1
end
cb = ContinuousCallback(cond, affect!)
@info "Event Repeat Test 1"
sol = solve(prob, Tsit5(); adaptive = false, callback = cb, dt = dt, save_everystep = false)
@test c[] == 1

function condition_v(out, u, t, integrator)
    out[1] = u - exp(50eps(1.0))
    out[2] = u - exp(60eps(1.0))
    out[3] = u - exp(70eps(1.0))
    out[4] = u - exp(80eps(1.0))
    return out[5] = u - exp(90eps(1.0))
end

c1 = Ref(0)
c2 = Ref(0)
c3 = Ref(0)
c4 = Ref(0)
c5 = Ref(0)

function affect_v!(integrator, idx)
    return if idx == 1
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

cb = VectorContinuousCallback(condition_v, affect_v!, 5)
@info "Event Repeat Test 2"
sol = solve(prob, Tsit5(); adaptive = false, callback = cb, dt = dt, save_everystep = false)

@test c1[] == 1
@test c2[] == 1
@test c3[] == 1
@test c4[] == 1
@test c5[] == 1
