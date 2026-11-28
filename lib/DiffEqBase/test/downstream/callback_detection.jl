using OrdinaryDiffEq, SciMLBase
import DiffEqBase, ReverseDiff
# https://github.com/SciML/DiffEqBase.jl/issues/1231
@testset "Successive different callbacks in same integration step" begin
    cb = ContinuousCallback(
        (u, t, integrator) -> t - 0.0,
        (integrator) -> push!(record, 0);
        abstol = 0.0
    )

    vcb = VectorContinuousCallback(
        (out, u, t, integrator) -> out .= (t - 1.0e-8, t - 2.0e-8, t - 2.0e-7),
        (integrator, events) -> begin
            for (idx, dir) in enumerate(events)
                iszero(dir) && continue
                push!(record, idx)
            end
        end,
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
    affect1 = (integ) -> push!(record, 1)
    affect2 = (integ) -> push!(record, 2)
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

@testset "Discontinuity detection is opt-in" begin
    prob = ODEProblem((u, p, t) -> u, 1.0, (0.0, 2.0))

    function solve_with_callback(maybe_discontinuity)
        callback = ContinuousCallback(
            (u, t, integrator) -> t - 0.5,
            integrator -> nothing;
            maybe_discontinuity
        )
        return solve(prob, Tsit5(); callback)
    end

    sol_default = solve_with_callback(false)
    sol_discontinuity_candidate = solve_with_callback(true)
    @test sol_default.t == sol_discontinuity_candidate.t
end

@testset "Successive same event detection" begin
    @testset "affect_integrator: affect_integrator" for affect_integrator in [false, true]
        @testset "tdir: $tdir" for tdir in [1, -1]
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
                affect!; abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = cb)
            @test record == [0.1, 0.4, 0.8]

            vcb = VectorContinuousCallback(
                (out, u, t, integrator) -> out .= (poly(tdir * t), poly(tdir * t - 0.1)),
                affect!, 2; abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = vcb)
            @test record == [0.1, 0.2, 0.4, 0.5, 0.8, 0.9]

            # Quadratic roots (cannot step on exact root)

            cb = ContinuousCallback(
                (u, t, integrator) -> poly(t^2),
                affect!; abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = cb)
            @test record ≈ sqrt.([0.1, 0.4, 0.8])

            vcb = VectorContinuousCallback(
                (out, u, t, integrator) -> out .= (poly(t^2), poly(t^2 - 0.1)),
                affect!, 2; abstol
            )

            record = []
            sol = solve(prob, Tsit5(), callback = vcb)
            @test record ≈ sqrt.([0.1, 0.2, 0.4, 0.5, 0.8, 0.9])
        end
    end
end

@testset "ReverseDiff at a sharp continuous callback" begin
    function drift!(du, u, p, t)
        du[1] = p[2]
        return nothing
    end
    condition(u, t, integrator) = atan((only(u) - integrator.p[1]) / 1.0e-6)
    vector_condition!(out, u, t, integrator) = (out[1] = condition(u, t, integrator))
    function reset!(integrator, args...)
        integrator.u[1] += 2
        return terminate!(integrator)
    end
    p = [0.43, 1.0]
    for rootfind in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind),
            callback in (
                ContinuousCallback(condition, reset!; rootfind),
                VectorContinuousCallback(vector_condition!, reset!, 1; rootfind),
            )
        function event_solve(x)
            u0 = [convert(eltype(x), 0.07)]
            tspan = convert.(eltype(x), (0.0, 1.0))
            prob = ODEProblem{true, SciMLBase.FullSpecialize}(drift!, u0, tspan, x)
            return solve(
                prob, Tsit5(); callback, save_everystep = false,
                sensealg = DiffEqBase.SensitivityADPassThrough()
            )
        end
        # The reset gives u = p[1] + 2 at t = (p[1] - 0.07) / p[2].
        @test only(last(event_solve(p).u)) ≈ 2.43
        @test ReverseDiff.gradient(x -> only(last(event_solve(x).u)), p) ≈ [1.0, 0.0] atol = 1.0e-6
        @test ReverseDiff.gradient(x -> last(event_solve(x).t), p) ≈ [1.0, -0.36] atol = 1.0e-6
    end
end

@testset "ReverseDiff replay of public continuous callbacks" begin
    function stopped_solve(x, condition, span, initial, dt, rootfind, vector)
        condition!(out, u, t, integrator) = (out[1] = condition(u, t, integrator))
        function stop!(integrator, args...)
            integrator.u += 2
            return terminate!(integrator)
        end
        callback = if vector
            VectorContinuousCallback(condition!, stop!, 1; rootfind, interp_points = 0)
        else
            ContinuousCallback(condition, stop!; rootfind, interp_points = 0)
        end
        prob = ODEProblem((u, p, t) -> one(u), initial, span, (x[1], x[2]))
        return solve(prob, Tsit5(); callback, dt, adaptive = false, sensealg = DiffEqBase.SensitivityADPassThrough())
    end

    left = SciMLBase.LeftRootFind
    right = SciMLBase.RightRootFind
    function cubic(u, t, integrator)
        q, k = integrator.p
        s = t - q
        return s * (1 + k * s^2)
    end
    function confined(u, t, integrator)
        q, w = integrator.p
        s = t - q
        return s * (1 - log1p(-(s / (2w))^2))
    end
    cases = (
        (cubic, x -> 0.25, [0.25, 1.0], ([0.3125, 1.0], [0.375, 1.0]), false, -1, right),
        (cubic, x -> 0.25, [0.5, 1.0], ([0.625, 1.0e14],), true, 1, left),
        (confined, last, [0.5, 0.25], ([0.5, 2.0^-21],), false, 1, right),
        (confined, last, [0.5, 0.25], ([0.5, 2.0^-21],), true, -1, left),
    )
    function replay_gradient(f, recorded, inputs, exact_value, gradient; value_atol, grad_atol)
        @test f(first(inputs)) ≈ exact_value(first(inputs)) atol = value_atol
        tape = ReverseDiff.GradientTape(f, recorded)
        compiled = ReverseDiff.compile(tape)
        for x in inputs
            @test ReverseDiff.gradient(f, x) ≈ gradient(x) atol = grad_atol
            for trace in (tape, compiled)
                @test ReverseDiff.gradient!(trace, x) ≈ gradient(x) atol = grad_atol
            end
        end
    end
    for (condition, width, recorded, inputs, vector, dir, side) in cases
        function event_squared(x)
            w = width(x)
            q = first(x)
            sol = stopped_solve(x, condition, (q - dir * w, q + dir * w), zero(q), dir * 2w, side, vector)
            return last(sol.t)^2
        end
        replay_gradient(
            event_squared, recorded, inputs, x -> first(x)^2, x -> [2first(x), 0.0];
            value_atol = 1.0e-10, grad_atol = 5.0e-6
        )
    end

    function localized(u, t, integrator)
        q, w = integrator.p
        return (u - q) * (1 + exp(-((u - 0.5) / w)^2))
    end
    for (vector, dir, side) in ((false, 1, left), (true, -1, right)), saved_state in (false, true)
        function event_output(x)
            z = zero(first(x))
            span = dir > 0 ? (z, z + 1) : (z + 1, z)
            sol = stopped_solve(x, localized, span, first(span), dir * 0.75, side, vector)
            return last(saved_state ? sol.u : sol.t)
        end
        value = saved_state ? 2.5 : 0.5
        @test event_output([0.5, 1.0e-8]) ≈ value atol = 2.0e-8
        @test ReverseDiff.gradient(event_output, [0.5, 1.0e-8]) ≈ [1.0, 0.0] atol = 2.0e-6
        replay_gradient(
            event_output, [0.5, 0.01], ([0.5, 1.0e-6],), x -> value, x -> [1.0, 0.0];
            value_atol = 2.0e-8, grad_atol = 2.0e-6
        )
    end
end
