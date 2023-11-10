using OrdinaryDiffEq, Test

# Setup a simple ODE problem with several callbacks (to test LLVM code gen)
# We will manually trigger the first callback and check its allocations.
function f!(du, u, p, t)
    du .= .-u
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
    integrator.p[1] += 1
end

cbs = CallbackSet(ContinuousCallback(cond_1, cb_affect!),
                  ContinuousCallback(cond_2, cb_affect!),
                  ContinuousCallback(cond_3, cb_affect!),
                  ContinuousCallback(cond_4, cb_affect!),
                  ContinuousCallback(cond_5, cb_affect!),
                  ContinuousCallback(cond_6, cb_affect!),
                  ContinuousCallback(cond_7, cb_affect!),
                  ContinuousCallback(cond_8, cb_affect!),
                  ContinuousCallback(cond_9, cb_affect!))

integrator = init(ODEProblem(f!, [0.8, 1.0], (0.0, 100.0), [0, 0]), Tsit5(), callback = cbs,
                  save_on = false);
# Force a callback event to occur so we can call handle_callbacks! directly.
# Step to a point where u[1] is still > 0.5, so we can force it below 0.5 and
# call handle callbacks
step!(integrator, 0.1, true)

if VERSION >= v"1.7"
    function handle_allocs(integrator)
        integrator.u[1] = 0.4
        @allocations OrdinaryDiffEq.handle_callbacks!(integrator)
    end
    handle_allocs(integrator);
    @test handle_allocs(integrator) == 0
end
