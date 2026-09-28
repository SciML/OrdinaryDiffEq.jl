using OrdinaryDiffEqSDIRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF
using OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner, SciMLBase, Test

@testset "Callback interpolation on a tstop-clamped step" begin
    for alg in (ImplicitEuler(), TRBDF2(), Rodas5P(), QBDF2(), QBDF(), FBDF(), Tsit5(), Vern9()),
            inplace in (false, true), direction in (1.0, -1.0), rewind in (false, true)
        @testset "$(nameof(typeof(alg))) inplace=$inplace direction=$direction rewind=$rewind" begin
            done = Ref(false)
            midpoint = direction / 2
            f = inplace ? ((du, u, p, t) -> fill!(du, 1)) : ((u, p, t) -> one.(u))
            prob = ODEProblem(f, [1.0], (0.0, 2direction))
            callback = DiscreteCallback(
                (u, t, integrator) -> !done[] && t == direction,
                function (integrator)
                    @test integrator.tprev == 0
                    @test integrator(midpoint) ≈ [1 + midpoint]
                    @test integrator(midpoint; idxs = 1) ≈ 1 + midpoint
                    @test integrator(midpoint, Val{1}) ≈ [1.0]
                    out = similar(integrator.u)
                    integrator(out, midpoint)
                    @test out ≈ [1 + midpoint]
                    integrator(out, midpoint, Val{1})
                    @test out ≈ [1.0]
                    proposal = get_proposed_dt(integrator)
                    if rewind
                        change_t_via_interpolation!(integrator, midpoint)
                        @test integrator.t == midpoint
                        @test integrator.u ≈ [1 + midpoint]
                        @test get_proposed_dt(integrator) == proposal
                    end
                    done[] = true
                    return nothing
                end;
                save_positions = (false, false)
            )
            # BDF startup estimates can reject this analytically exact first step.
            startup = alg isa Union{QNDF, FBDF} ? (dtmin = 1.0, force_dtmin = true) : (;)
            sol = solve(
                prob, alg; dt = 5direction, tstops = [direction], callback,
                save_everystep = false, abstol = 1.0e-8, reltol = 1.0e-8, startup...
            )
            @test done[]
            @test SciMLBase.successful_retcode(sol)
            @test sol.u[end] ≈ [1 + 2direction]
        end
    end
end

# On a first step clamped onto a tstop, saveat values and ContinuousCallback
# event times must use the accepted interval (exact solution u = 1 + t).
@testset "saveat and ContinuousCallback on a tstop-clamped first step" begin
    for alg in (Rodas5P(), ImplicitEuler(), TRBDF2(), Tsit5(), Vern9())
        @testset "$(nameof(typeof(alg)))" begin
            prob = ODEProblem((u, p, t) -> one.(u), [1.0], (0.0, 2.0))
            sol = solve(
                prob, alg; dt = 5.0, tstops = [1.0], saveat = [0.5, 1.5],
                abstol = 1.0e-8, reltol = 1.0e-8
            )
            @test SciMLBase.successful_retcode(sol)
            @test sol.t == [0.5, 1.5]
            @test sol.u[1] ≈ [1.5]
            @test sol.u[2] ≈ [2.5]

            event_t = Ref(NaN)
            cb = ContinuousCallback(
                (u, t, integrator) -> u[1] - 1.5,
                function (integrator)
                    event_t[] = integrator.t
                    return nothing
                end;
                save_positions = (false, false)
            )
            sol_cb = solve(
                prob, alg; dt = 5.0, tstops = [1.0], callback = cb,
                save_everystep = false, abstol = 1.0e-8, reltol = 1.0e-8
            )
            @test SciMLBase.successful_retcode(sol_cb)
            @test event_t[] ≈ 0.5
        end
    end
end
