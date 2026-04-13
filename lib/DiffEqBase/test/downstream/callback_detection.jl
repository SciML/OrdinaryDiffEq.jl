using OrdinaryDiffEq
# https://github.com/SciML/DiffEqBase.jl/issues/1231
@testset "Successive different callbacks in same integration step" begin
    cb = ContinuousCallback(
        (u, t, integrator) -> t - 0.0,
        (integrator) -> push!(record, 0);
        abstol = 0.0
    )

    vcb = VectorContinuousCallback(
        (out, u, t, integrator) -> out .= (t - 1.0e-8, t - 2.0e-8, t - 2.0e-7),
        (integrator, event_index) -> push!(record, event_index),
        3;
        abstol = 0.0
    )

    f(u, p, t) = 1.0
    u0 = 0.0

    # Forward propagation with successive events
    record = []
    tspan = (-1.0, 1.0)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Tsit5(), dt = 2.0, callback = CallbackSet(cb, vcb))
    @test record == [0, 1, 2, 3]

    # Backward propagation with successive events
    record = []
    tspan = (1.0, -1.0)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Tsit5(), dt = 2.0, callback = CallbackSet(cb, vcb))
    @test record == [3, 2, 1, 0]

    # https://github.com/SciML/DiffEqBase.jl/issues/1268
    # Test that callback order does not affect event detection
    function rhs!(du, u, p, t)
        du[1] = -3sin(3t)
        du[2] = -sin(t)
    end
    affect1 = (integ,) -> push!(record, 1)
    affect2 = (integ,) -> push!(record, 2)
    cb1 = ContinuousCallback((u, i, integ) -> u[2], affect1, affect1; rootfind = SciMLBase.LeftRootFind, abstol = 0.0)
    cb2 = ContinuousCallback((u, i, integ) -> u[1], affect2, affect2; rootfind = SciMLBase.RightRootFind, abstol = 0.0)

    record = []
    prob = ODEProblem(rhs!, ones(2), (0.0, 2π); callback = CallbackSet(cb2, cb1))
    sol = solve(prob, Tsit5())
    @test record == [2, 2, 1, 2, 2, 1, 2, 2]

    record = []
    prob = ODEProblem(rhs!, ones(2), (0.0, 2π); callback = CallbackSet(cb1, cb2))
    sol = solve(prob, Tsit5())
    @test record == [2, 2, 1, 2, 2, 1, 2, 2]
end

@testset "Successive same event detection" begin
    @testset for affect_integrator in [false, true]
        @testset for tdir in [1, -1]
            poly(t) = (t - 0.1) * (t - 0.4) * (t - 0.8)
            function affect!(integrator, index = 1)
                push!(record, tdir * integrator.t)
                if affect_integrator
                    # nudge t backward to see if integrator avoids repeat detection
                    integrator.t = integrator.t - tdir * 1.0e-14
                end
            end
            abstol = affect_integrator ? 1.0e-14 : 0.0

            f(u, p, t) = 1.0
            u0 = 0.0
            tspan = tdir .* (0.0, 1.0)
            prob = ODEProblem(f, u0, tspan; dt = 0.25, maxiters = 100)

            # Linear roots (can step on exact root)

            cb = ContinuousCallback(
                (u, t, integrator) -> poly(tdir * t),
                affect!; abstol = abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = cb)
            @test record == [0.1, 0.4, 0.8]

            vcb = VectorContinuousCallback(
                (out, u, t, integrator) -> out .= (poly(tdir * t), poly(tdir * t - 0.1)),
                affect!, 2; abstol = abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = vcb)
            @test record == [0.1, 0.2, 0.4, 0.5, 0.8, 0.9]

            # Quadratic roots (cannot step on exact root)

            cb = ContinuousCallback(
                (u, t, integrator) -> poly(t^2),
                affect!; abstol = abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = cb)
            @test record ≈ sqrt.([0.1, 0.4, 0.8])

            vcb = VectorContinuousCallback(
                (out, u, t, integrator) -> out .= (poly(t^2), poly(t^2 - 0.1)),
                affect!, 2; abstol = abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = vcb)
            @test record ≈ sqrt.([0.1, 0.2, 0.4, 0.5, 0.8, 0.9])
        end
    end
end
