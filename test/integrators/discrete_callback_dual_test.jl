using OrdinaryDiffEq, Test, ForwardDiff

# Test problem
# u' = p(t)u; p(t) = p if t <= tstop, 0 otherwise
# Target: du(1.0)/dtstop = exp(tstop)u0

u0 = 1.0
tspan = (0.0, 1.0)
p = 1.0

function stopping_cb(tstop)
    condition = (u, t, integrator) -> t == tstop
    affect! = integrator -> (println("Stopped!"); integrator.p = zero(integrator.p))
    DiscreteCallback(condition, affect!)
end

function test_fun(tstop)
    DualT = typeof(tstop)
    prob = ODEProblem((u, p, t) -> p * u, DualT(u0), DualT.(tspan), DualT(p))
    sol = solve(prob, Tsit5(), callback = stopping_cb(tstop), tstops = [tstop])
    sol(1.0)
end

@test ForwardDiff.derivative(test_fun, 0.5) ≈ exp(0.5) * u0 # Analytical solution: exp(tstop)*u0

function test_fun(tstop)
    DualT = typeof(tstop)
    prob = ODEProblem((u, p, t) -> p * u, DualT(u0), DualT.(tspan), DualT(p))
    sol = solve(prob, Tsit5(), callback = stopping_cb(tstop), tstops = [tstop],
        adaptive = false, dt = 0.01)
    sol(1.0)
end

@test ForwardDiff.derivative(test_fun, 0.5) ≈ exp(0.5) * u0 # Analytical solution: exp(tstop)*u0
