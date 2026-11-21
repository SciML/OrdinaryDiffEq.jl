using OrdinaryDiffEqBDF, Test
using OrdinaryDiffEqCore: set_EEst!

@testset "MOOSE234 variable-step order-four update" begin
    for iip in (false, true), mesh in ((0, -1, -2, -3), (0, -2, -3, -4)), h in (0.125, 0.0625)
        f = iip ? ((du, u, p, t) -> (du[1] = 4t^3)) : ((u, p, t) -> 4t^3)
        prob = ODEProblem(f, iip ? [0.0] : 0.0, (0.0, 1.0))
        integ = init(prob, MOOSE234(); dt = h, abstol = 1.0e-10, reltol = 1.0e-10)
        cache = integ.cache
        cache.order = 4
        cache.iters_from_event = 3
        cache.consfailcnt = 1
        for (i, offset) in enumerate(mesh)
            cache.ts[i] = h * offset
            value = cache.ts[i]^4
            if iip
                cache.u_history[i][1] = value
            else
                cache.u_history[i] = value
            end
        end
        OrdinaryDiffEqBDF.perform_step!(integ, cache)
        @test !integ.force_stepfail
        @test (iip ? only(integ.u) : integ.u) ≈ h^4 rtol = 1.0e-10
    end
end

@testset "MOOSE234 startup interpolation and domain retry" begin
    for iip in (false, true)
        f = iip ? ((du, u, p, t) -> fill!(du, 1.0)) : ((u, p, t) -> 1.0)
        prob = ODEProblem(f, iip ? [1.0] : 1.0, (0.0, 0.8))
        scalar = iip ? only : identity
        callback = ContinuousCallback((u, t, integ) -> t - 0.25, terminate!)
        sol = solve(prob, MOOSE234(); callback, dt = 0.5, abstol = 1.0, reltol = 0.1)
        @test sol.retcode == ReturnCode.Terminated
        @test sol.t[end] ≈ 0.25
        @test scalar(sol(0.125)) ≈ 1.125

        trials = Ref(0)
        isoutofdomain(u, p, t) = (trials[] += 1; trials[] == 3)
        sol = solve(
            prob, MOOSE234(); isoutofdomain, dt = 0.1, dtmax = 0.1,
            abstol = 1.0, reltol = 0.1
        )
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject == 1
        @test scalar(sol.u[end]) ≈ 1.8
    end
end

@testset "MOOSE234 controller time direction" begin
    for iip in (false, true)
        f = iip ? ((du, u, p, t) -> fill!(du, 1.0)) : ((u, p, t) -> 1.0)
        for direction in (1, -1)
            prob = ODEProblem(f, iip ? [1.0] : 1.0, (0.0, Float64(direction)))
            integ = init(prob, MOOSE234(); dt = direction * 0.1)
            cache = integ.cache
            cache.iters_from_event = 3
            cache.qwait = 0
            cache.terkm1 = 0.1
            cache.terkp1 = 0.5
            @test first(OrdinaryDiffEqBDF.choose_order_moose!(integ, cache)) == 2

            cache.order = 3
            set_EEst!(integ, 2.0)
            OrdinaryDiffEqBDF.step_reject_controller!(integ, integ.alg)
            @test cache.order == 2
            @test integ.dt ≈ direction * 0.1
        end
    end
end
