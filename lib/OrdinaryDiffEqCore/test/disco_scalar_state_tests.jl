using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, Test

noop!(integrator) = nothing
disco_controller(alg) = OrdinaryDiffEqCore.PIController(alg; discontinuity_detection = true)

@testset "Discontinuity detection with non-mutable states" begin
    @testset "find_discontinuity with a scalar out-of-place state" begin
        prob = ODEProblem((u, p, t) -> 1.0, 0.0, (0.0, 2.0))
        cb = ContinuousCallback((u, t, integrator) -> t - 0.3, noop!)
        integrator = init(
            prob, Tsit5(); callback = cb, dt = 0.1, adaptive = false,
            controller = disco_controller(Tsit5())
        )
        integrator.dt = 1.0
        integrator.u = 1.0
        θ = OrdinaryDiffEqCore.find_discontinuity(integrator)
        @test θ ≈ 0.3 atol = 1.0e-8
        @test integrator.is_disco_step
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
