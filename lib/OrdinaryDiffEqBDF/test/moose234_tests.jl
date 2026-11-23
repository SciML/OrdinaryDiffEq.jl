using OrdinaryDiffEqBDF, Test
using OrdinaryDiffEqCore: set_EEst!
using LinearAlgebra: Diagonal, Hermitian, Symmetric

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

@testset "BDF zero-valued dense sample" begin
    zero_time = (1 + cospi(3 / 4)) / 2
    for alg in (MOOSE234(), FBDF()), iip in (false, true)
        f = iip ? ((du, u, p, t) -> fill!(du, 1.0)) : ((u, p, t) -> 1.0)
        prob = ODEProblem(f, iip ? [-zero_time] : -zero_time, (0.0, 1.0))
        scalar = iip ? only : identity
        opts = (; dt = 1.0, abstol = 1.0, reltol = 1.0e-6)
        sol = solve(prob, alg; opts...)
        @test scalar(sol(0.5)) ≈ 0.5 - zero_time
        @test scalar(sol(0.5, Val{1})) ≈ 1.0
        saved = solve(prob, alg; opts..., saveat = [0.5])
        @test scalar(only(saved.u)) ≈ 0.5 - zero_time
        callback = ContinuousCallback((u, t, integ) -> scalar(u), terminate!)
        event = solve(prob, alg; opts..., callback)
        @test event.retcode == ReturnCode.Terminated
        @test event.t[end] ≈ zero_time
        @test scalar(event(zero_time / 2)) ≈ -zero_time / 2
    end
end

@testset "MOOSE234 singular mass matrix" begin
    mass_matrix = [1.0 0.0; 0.0 0.0]
    for wrap in (Diagonal, Symmetric, Hermitian), iip in (false, true)
        f = iip ? ((du, u, p, t) -> (du .= (1.0, u[2] - exp(u[1])))) :
            ((u, p, t) -> [1.0, u[2] - exp(u[1])])
        prob = ODEProblem(ODEFunction(f; mass_matrix = wrap(mass_matrix)), [0.0, 1.0], (0.0, 2.0))
        @test_throws ArgumentError("MOOSE234 requires a nonsingular mass matrix.") init(
            prob, MOOSE234(); abstol = 1.0e-6, reltol = 1.0e-6
        )
    end
end

@testset "MOOSE234 invertible indefinite mass matrix" begin
    mass_matrix = [0.0 1.0; 1.0 0.0]
    for wrap in (Symmetric, Hermitian), iip in (false, true)
        @testset "$wrap, in-place = $iip" begin
            f = iip ? ((du, u, p, t) -> (du .= -u)) : ((u, p, t) -> -u)
            u0 = [1.0, 1.0]
            prob = ODEProblem(ODEFunction(f; mass_matrix = wrap(mass_matrix)), u0, (0.0, 1.0))
            sol = solve(prob, MOOSE234(); dt = 0.001, abstol = 1.0e-8, reltol = 1.0e-8)
            @test sol.retcode == ReturnCode.Success
            @test sol.u[end] ≈ exp(-1) .* u0 atol = 1.0e-8 rtol = 0
        end
    end
end

@testset "MOOSE234 step growth and safety factor" begin
    for iip in (false, true)
        f = iip ? ((du, u, p, t) -> (du .= -u)) : ((u, p, t) -> -u)
        prob = ODEProblem(f, iip ? [1.0] : 1.0, (0.0, 1.0))
        opts = (; dt = 1.0e-4, abstol = 1.0e-6, reltol = 1.0e-6)
        controller = OrdinaryDiffEqBDF.BDFController(qmax = 1.1, qmax_first_step = 1.1)
        sol = solve(prob, MOOSE234(); opts..., controller)
        @test sol.retcode == ReturnCode.Success
        steps = diff(sol.t)
        @test maximum(steps[2:end] ./ steps[1:(end - 1)]) ≈ 1.1 rtol = 1.0e-12

        cautious = solve(prob, MOOSE234(); opts..., controller = OrdinaryDiffEqBDF.BDFController(gamma = 0.25))
        aggressive = solve(prob, MOOSE234(); opts..., controller = OrdinaryDiffEqBDF.BDFController(gamma = 0.9))
        @test cautious.retcode == aggressive.retcode == ReturnCode.Success
        @test cautious.t[5] < aggressive.t[5]
    end
end
