using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, Test

noop!(integrator) = nothing
disco_controller(alg) = OrdinaryDiffEqCore.PIController(alg; discontinuity_detection = true)

@testset "Discontinuity detection with non-mutable states" begin
    @testset "find_discontinuity with a scalar out-of-place state ($T, $spec)" for T in (Float64, Float32),
            spec in (SciMLBase.FullSpecialize, SciMLBase.AutoSpecialize)
        prob = ODEProblem{false, spec}((u, p, t) -> 1.0, 0.0, (zero(T), T(2)))
        cb = ContinuousCallback((u, t, integrator) -> t - T(0.3), noop!)
        integrator = init(
            prob, Tsit5(); callback = cb, dt = T(0.1), adaptive = false,
            controller = disco_controller(Tsit5())
        )
        integrator.dt = one(T)
        integrator.u = 1.0
        θ = OrdinaryDiffEqCore.find_discontinuity(integrator)
        @test θ ≈ T(0.3) atol = (T === Float64 ? 1.0e-8 : 4eps(0.3f0))
        @test integrator.is_disco_step
        probs = OrdinaryDiffEqCore.get_disco_probs(integrator.controller_cache)
        @test isconcretetype(eltype(probs))
        disco = only(probs)
        @test disco.f.f isa (spec === SciMLBase.FullSpecialize ? OrdinaryDiffEqCore.zero_func_struct : OrdinaryDiffEqCore.FunctionWrapper)
        @test disco.f(T(0.1), disco.p) === Float64(T(0.1)) - Float64(T(0.3))
    end

    # `maybe_discontinuity` does not change a callback's type, so a scan with every root skipped
    # runs the same code without the root solves. Each solve must then cost the same in a mixed
    # set as for one scalar callback (nothing, once the bracketing solver does not splat the
    # bracket). The roots decrease along each set (vector 3/8 and 1/4, then 0.2, then 0.1), so
    # every bracket encloses its root and each crossing is solved once.
    scalar(md) = ContinuousCallback((u, t, integrator) -> t - 0.2, noop!; maybe_discontinuity = md)
    scalar2(md) = ContinuousCallback((u, t, integrator) -> t - 0.1, noop!; maybe_discontinuity = md)
    vector(md) = VectorContinuousCallback(
        (out, u, t, integrator) -> (out[1] = t - 3 / 8; out[2] = t - 1 / 4; nothing),
        (integrator, idx) -> nothing, 2; maybe_discontinuity = md
    )
    function scan(cb)
        integrator = init(
            ODEProblem{true, SciMLBase.FullSpecialize}((du, u, p, t) -> (du[1] = 1.0; nothing), [0.0], (0.0, 2.0)),
            Tsit5(); callback = cb, dt = 0.5, adaptive = false
        )
        integrator.dt = 1.0
        integrator.u .= integrator.uprev .+ 1.0
        θ = OrdinaryDiffEqCore.find_discontinuity(integrator)
        return θ, @allocated OrdinaryDiffEqCore.find_discontinuity(integrator)
    end
    function solve_allocations(cbs)
        θ, solved = scan(cbs(true))
        θ_skipped, skipped = scan(cbs(false))
        @test θ_skipped == 1
        return θ, solved - skipped
    end
    @testset "Discontinuity scan allocations under FullSpecialize: $name" for (name, cbs, root, nsolves) in (
            ("one vector callback", vector, 0.25, 2),
            ("vector-backed callback set", md -> CallbackSet([vector(md)], DiscreteCallback[]), 0.25, 2),
            ("vector and an ineligible scalar callback", md -> CallbackSet(vector(md), scalar(false)), 0.25, 2),
            ("vector and two scalar callbacks", md -> CallbackSet(vector(md), scalar(md), scalar2(md)), 0.1, 4),
        )
        θ_scalar, per_solve = solve_allocations(scalar)
        @test θ_scalar ≈ 0.2 atol = 1.0e-8
        θ, allocations = solve_allocations(cbs)
        @test θ ≈ root atol = 1.0e-8
        @test allocations == nsolves * per_solve
    end

    # A callback vector whose element type gives no concrete problem type keeps the wrapper, so the
    # storage stays concrete rather than abstract.
    @testset "Abstract callback vector under FullSpecialize: $name" for (name, cbs, root) in (
            ("one vector callback", SciMLBase.AbstractContinuousCallback[vector(true)], 0.25),
            ("vector and scalar callbacks", SciMLBase.AbstractContinuousCallback[vector(true), scalar(true)], 0.2),
        )
        integrator = init(
            ODEProblem{true, SciMLBase.FullSpecialize}((du, u, p, t) -> (du[1] = 1.0; nothing), [0.0], (0.0, 2.0)),
            Tsit5(); callback = CallbackSet(cbs, DiscreteCallback[]), dt = 0.5, adaptive = false
        )
        integrator.dt = 1.0
        integrator.u .= integrator.uprev .+ 1.0
        @test OrdinaryDiffEqCore.find_discontinuity(integrator) ≈ root atol = 1.0e-8
        probs = OrdinaryDiffEqCore.get_disco_probs(integrator.controller_cache)
        @test isconcretetype(eltype(probs))
        @test all(prob -> prob.f.f isa OrdinaryDiffEqCore.FunctionWrapper, probs)
    end

    # u' = 2u below u = 1 and u' = 5 - 3u above it, from u(0) = 0.1: the crossing is at
    # τ = log(10)/2 and u(t) = 5/3 - (2/3)exp(-3(t - τ)) afterwards.
    kink(u) = u < 1 ? 2u : 5 - 3u
    τ_kink = log(10) / 2
    u_kink(t) = t < τ_kink ? 0.1 * exp(2t) : 5 / 3 - (2 / 3) * exp(-3(t - τ_kink))

    @testset "Discontinuity detection solve, state type $(label)" for (label, prob, getu) in (
            ("Float64", ODEProblem((u, p, t) -> kink(u), 0.1, (0.0, 1.5)), identity),
            ("Vector, out-of-place", ODEProblem((u, p, t) -> [kink(u[1])], [0.1], (0.0, 1.5)), first),
        )
        zero_func_calls = Ref(0)
        function condition(u, t, integrator)
            integrator isa OrdinaryDiffEqCore.zero_func_struct && (zero_func_calls[] += 1)
            return getu(u) - 1
        end
        sol = solve(
            prob, Tsit5(); callback = ContinuousCallback(condition, noop!), dt = 1.4,
            reltol = 1.0e-8, abstol = 1.0e-10, controller = disco_controller(Tsit5())
        )
        @test sol.retcode == ReturnCode.Success
        @test zero_func_calls[] > 0
        @test getu(sol.u[end]) ≈ u_kink(1.5) rtol = 1.0e-6
        @test getu(sol(1.0)) ≈ u_kink(1.0) rtol = 1.0e-6
    end

    # u' = 1 before t = 0.3 and u' = -1 after it, from u(0) = 0: u(t) = 0.6 - t for t ≥ 0.3.
    @testset "Time-only discontinuity, scalar out-of-place state" begin
        prob = ODEProblem((u, p, t) -> t < 0.3 ? 1.0 : -1.0, 0.0, (0.0, 1.0))
        zero_func_calls = Ref(0)
        function condition(u, t, integrator)
            integrator isa OrdinaryDiffEqCore.zero_func_struct && (zero_func_calls[] += 1)
            return t - 0.3
        end
        sol = solve(
            prob, Tsit5(); callback = ContinuousCallback(condition, noop!), dt = 0.9,
            reltol = 1.0e-8, abstol = 1.0e-10, controller = disco_controller(Tsit5())
        )
        @test sol.retcode == ReturnCode.Success
        @test zero_func_calls[] > 0
        @test sol.u[end] ≈ -0.4 atol = 1.0e-6
        @test sol(0.3) ≈ 0.3 atol = 1.0e-6
    end
end
