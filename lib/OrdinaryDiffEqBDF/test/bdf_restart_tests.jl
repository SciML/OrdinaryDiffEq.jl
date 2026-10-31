using OrdinaryDiffEqBDF, SciMLBase, Test

@testset "BDF history after a discontinuity" begin
    for (name, alg) in ("QNDF" => QNDF(), "QNDF2" => QNDF2(), "QBDF" => QBDF(), "QBDF2" => QBDF2()),
            inplace in (false, true),
            direction in (1.0, -1.0), adaptive in (false, true), modification in (:state, :rewind)
        @testset "$name inplace=$inplace direction=$direction adaptive=$adaptive modification=$modification" begin
            f = inplace ? ((du, u, p, t) -> (du .= p[1] .* u .+ p[2])) :
                ((u, p, t) -> p[1] .* u .+ p[2])
            prob = ODEProblem(f, [1.0], (0.0, 3direction), [1.0, 0.0])
            integrator = init(prob, alg; dt = direction / 8, dtmax = 1 / 8, adaptive, save_everystep = false, abstol = 1.0e-10, reltol = 1.0e-10)
            while direction * integrator.t < 1
                step!(integrator)
            end
            if modification == :rewind
                change_t_via_interpolation!(integrator, (integrator.tprev + integrator.t) / 2)
            end
            event_time = integrator.t
            slope = startswith(name, "QBDF") ? 1.0 : 0.0
            integrator.u .= 4
            integrator.p .= (0, slope)
            derivative_discontinuity!(integrator, true)
            for _ in 1:8
                step!(integrator)
                @test integrator.u ≈ [4 + slope * (integrator.t - event_time)]
            end
        end
    end
end

# No-event QNDF2 with a user dt that triggers an early rejection must keep the
# attempt-based startup counter (`cnt` counts attempts since the last history
# reset; the first two use BDF1 coefficients).
@testset "QNDF2 no-event early rejection step count" begin
    function rober!(du, u, p, t)
        y₁, y₂, y₃ = u
        k₁, k₂, k₃ = p
        du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
        du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
        du[3] = k₂ * y₂^2
        return nothing
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e3), (0.04, 3.0e7, 1.0e4))
    sol = solve(prob, QNDF2(); dt = 1.0e-3, abstol = 1.0e-4, reltol = 1.0e-4)
    @test sol.retcode == ReturnCode.Success
    @test sol.stats.naccept == 54
    @test sol.stats.nreject == 1
end

# `reinit!` must not leave QNDF2 event anchors stale relative to reset `iter`,
# and must not wipe cold-start history on `initialize!` (see #4655).
@testset "QNDF2 reinit! preserves Success on Robertson" begin
    function rober!(du, u, p, t)
        y₁, y₂, y₃ = u
        k₁, k₂, k₃ = p
        du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
        du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
        du[3] = k₂ * y₂^2
        return nothing
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e3), (0.04, 3.0e7, 1.0e4))
    integ = init(prob, QNDF2(); dt = 1.0e-3, abstol = 1.0e-4, reltol = 1.0e-4)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success
    reinit!(integ)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success
end

# After a mid-solve event, `reinit!(; reinit_cache=false)` zeroes `iter` without
# calling `initialize!`, so `perform_step!` must re-anchor when
# `iter <= iter_at_event`. Stats are not reset by `reinit!`, so a correct
# re-solve doubles the first run's accept/reject counts (without the branch:
# Robertson QNDF2 goes ~130/5 → ~400/124).
@testset "QNDF2 reinit!(; reinit_cache=false) re-anchors after an event" begin
    function rober!(du, u, p, t)
        y₁, y₂, y₃ = u
        du[1] = -0.04y₁ + 1.0e4 * y₂ * y₃
        du[2] = 0.04y₁ - 3.0e7 * y₂^2 - 1.0e4 * y₂ * y₃
        du[3] = 3.0e7 * y₂^2
        return nothing
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 100.0))
    tev = 50.0
    cb = DiscreteCallback((u, t, integrator) -> t == tev, integrator -> (integrator.u .*= 1.0))
    integ = init(prob, QNDF2(); callback = cb, tstops = [tev], abstol = 1.0e-6, reltol = 1.0e-6)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success
    first_naccept = integ.stats.naccept
    first_nreject = integ.stats.nreject
    reinit!(integ; reinit_cache = false)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.Success
    @test integ.stats.naccept == 2 * first_naccept
    @test integ.stats.nreject == 2 * first_nreject
end
