using OrdinaryDiffEq, Test
using OrdinaryDiffEqCore
using OrdinaryDiffEqTsit5

# Setup a simple ODE problem with several callbacks (to test LLVM code gen)
# We will manually trigger the first callback and check its allocations.
function f!(du, u, p, t)
    return du .= .-u
end

cond_1(u, t, integrator) = u[1] - 0.5
cond_2(u, t, integrator) = u[2] + 0.5
cond_3(u, t, integrator) = u[2] + 0.6
cond_4(u, t, integrator) = u[2] + 0.7
cond_5(u, t, integrator) = u[2] + 0.8
cond_6(u, t, integrator) = u[2] + 0.9
cond_7(u, t, integrator) = u[2] + 0.1
cond_8(u, t, integrator) = u[2] + 0.11
cond_9(u, t, integrator) = u[2] + 0.12

function cb_affect!(integrator)
    return integrator.p[1] += 1
end

cbs = CallbackSet(
    ContinuousCallback(cond_1, cb_affect!),
    ContinuousCallback(cond_2, cb_affect!),
    ContinuousCallback(cond_3, cb_affect!),
    ContinuousCallback(cond_4, cb_affect!),
    ContinuousCallback(cond_5, cb_affect!),
    ContinuousCallback(cond_6, cb_affect!),
    ContinuousCallback(cond_7, cb_affect!),
    ContinuousCallback(cond_8, cb_affect!),
    ContinuousCallback(cond_9, cb_affect!)
)

integrator = init(
    ODEProblem{true, SciMLBase.FullSpecialize}(
        f!, [0.8, 1.0],
        (0.0, 100.0), [0, 0]
    ), Tsit5(), callback = cbs,
    save_on = false
);
# Force a callback event to occur so we can call handle_callbacks! directly.
# Step to a point where u[1] is still > 0.5, so we can force it below 0.5 and
# call handle callbacks
step!(integrator, 0.1, true)

function handle_allocs(integrator)
    integrator.u[1] = 0.4
    return @allocations OrdinaryDiffEqCore.handle_callbacks!(integrator)
end
handle_allocs(integrator)
@test_skip handle_allocs(integrator) == 0

@testset "AutoSpecialize erased callback step! allocations" begin
    noop!(integrator) = (derivative_discontinuity!(integrator, false); nothing)
    cond_never(u, t, integrator) = u[1] + 10.0
    f_alloc!(du, u, p, t) = (du .= .-u; nothing)
    u0 = [1.0, 0.0]
    tspan = (0.0, 1.0e6)
    prob = ODEProblem{true, SciMLBase.AutoSpecialize}(f_alloc!, u0, tspan)

    function measure_steps!(integ, n)
        return @allocated begin
            for _ in 1:n
                step!(integ)
            end
        end
    end
    function bytes_per_step(callback::CB; nwarm = 200, n = 1000) where {CB}
        integ = init(
            prob, Tsit5(); callback = callback,
            save_everystep = false, save_start = false, save_end = false, dense = false
        )
        for _ in 1:nwarm
            step!(integ)
        end
        measure_steps!(integ, 10)
        return measure_steps!(integ, n) / n
    end

    disc = DiscreteCallback((u, t, i) -> false, noop!; save_positions = (false, false))
    cont = ContinuousCallback(cond_never, noop!; save_positions = (false, false))

    @test bytes_per_step(disc) == 0
    @test bytes_per_step(cont) <= 64
end
