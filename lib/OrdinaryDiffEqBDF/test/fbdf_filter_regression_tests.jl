using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, SciMLBase, Test

function polynomial_history(k, h, spacing, iip; max_order = max(k + 1, 3), degree = k + 1)
    f(u, p, t) = degree * t^(degree - 1)
    f!(du, u, p, t) = (du[1] = f(u, p, t); nothing)
    prob = iip ? ODEProblem(f!, [0.0], (0.0, sign(h))) :
        ODEProblem(f, 0.0, (0.0, sign(h)))
    integ = init(
        prob, FBDF(time_filter = true, max_order = Val(max_order));
        dt = h, adaptive = false
    )
    cache = integ.cache
    cache.order = k
    cache.iters_from_event = k + 1
    cache.consfailcnt = 1
    integ.derivative_discontinuity = false
    for j in eachindex(cache.ts)
        cache.ts[j] = -(spacing isa Number ? (j - 1) * spacing : sum(spacing[1:(j - 1)]; init = 0.0)) * h
        value = cache.ts[j]^degree
        if iip
            cache.u_history[j] .= value
        else
            cache.u_history[j] = value
        end
    end
    return integ
end

@testset "FBDF filter regressions" begin
    @testset "Polynomial reproduction" begin
        for iip in (false, true), k in 1:4, spacing in (0.5, 1.0, 3.0, (0.5, 1.5, 0.75, 2.5, 1.25, 2.0)), h in (0.1, -0.1)
            integ = polynomial_history(k, h, spacing, iip)
            OrdinaryDiffEqBDF.perform_step!(integ, integ.cache)
            @test only(integ.u) ≈ h^(k + 1) rtol = 1.0e-9 atol = 1.0e-14
        end
    end
    @testset "Ineligible higher candidate" begin
        for iip in (false, true)
            integ = polynomial_history(2, 0.1, 1.0, iip)
            integ.opts.adaptive = true
            integ.cache.qwait = 3
            before = copy(integ.u)
            if iip
                moved = OrdinaryDiffEqBDF._fbdf_time_filter!(integ, integ.cache, 2)
                @test !moved
                @test integ.u == before
            else
                after, _ = OrdinaryDiffEqBDF._fbdf_time_filter(integ, integ.cache, before, 2)
                @test after == before
            end
        end
    end
    @testset "BDF3-Stab at max_order=3" begin
        for iip in (false, true)
            integ = polynomial_history(3, 0.1, 1.0, iip; max_order = 3)
            integ.opts.adaptive = true
            integ.cache.qwait = 0
            for j in eachindex(integ.cache.u_history)
                if iip
                    integ.cache.u_history[j] .= (j == 4 ? 1.0 : 0.0)
                else
                    integ.cache.u_history[j] = (j == 4 ? 1.0 : 0.0)
                end
            end
            if iip
                @test OrdinaryDiffEqBDF._fbdf_time_filter!(integ, integ.cache, 3)
            else
                OrdinaryDiffEqBDF._fbdf_time_filter(integ, integ.cache, integ.u, 3)
            end
            @test OrdinaryDiffEqCore.get_EEst(integ) == 0
        end
    end
end

@testset "FBDF candidate error and controller" begin
    for (k, max_order, qwait, errors, expected) in (
            (2, 5, 0, (0.5, 0.001, Inf), 3),
            (2, 5, 0, (0.001, 0.5, Inf), 2),
            (2, 5, 3, (0.5, 0.001, Inf), 2),
            (3, 5, 0, (0.5, 0.4, 0.001), 2),
            (3, 5, 0, (0.001, 0.5, 0.4), 3),
            (3, 5, 0, (0.5, 0.001, 0.4), 4),
            (3, 3, 0, (0.5, 0.001, 0.4), 2),
            (3, 5, 0, (0.5, NaN, 0.4), 2),
            (3, 5, 0, (0.0, 1.0e-30, 0.4), 3),
        )
        integ = polynomial_history(k, 0.1, 1.0, false; max_order)
        cache = integ.cache
        cache.qwait = qwait
        cache.nconsteps = 42
        integ.success_iter = 10
        selected = OrdinaryDiffEqBDF._fbdf_select_filter!(integ, cache, k, max_order, errors...)
        @test selected == expected
        @test OrdinaryDiffEqCore.get_current_alg_order(integ.alg, cache) == expected
        @test OrdinaryDiffEqCore.get_current_adaptive_order(integ.alg, cache) == expected
        err = errors[selected == k ? 1 : selected == k + 1 ? 2 : 3]
        @test OrdinaryDiffEqCore.get_EEst(integ) == err
        fill!(cache.stald.ssdat, 1)
        q = OrdinaryDiffEqBDF.stepsize_controller!(integ, cache, integ.alg)
        @test q ≈ clamp(err^(1 / (selected + 1)) * 1.2, 0.1, 5.0)
        @test all(iszero, cache.stald.ssdat)
        @test cache.order == k
        OrdinaryDiffEqBDF.step_accept_controller!(integ, cache, integ.alg, q)
        @test cache.order == max(k, expected)
        @test cache.nconsteps == (expected > k ? 1 : 43)
    end
end

@testset "FBDF filtered dense output and derivatives" begin
    for iip in (false, true), k in 1:4, spacing in (0.5, 3.0)
        h = 0.1
        integ = polynomial_history(k, h, spacing, iip)
        OrdinaryDiffEqBDF.perform_step!(integ, integ.cache)
        @test only(integ.fsallast) ≈ (k + 1) * h^k
        for θ in (0.2, 0.5, 0.8)
            y = OrdinaryDiffEqBDF._ode_interpolant(
                θ, h, integ.uprev, integ.u, integ.k, integ.cache, nothing, Val{0}, nothing
            )
            @test only(y) ≈ (θ * h)^(k + 1) rtol = 1.0e-8 atol = 1.0e-14
        end
    end
end

@testset "FBDF filtering during fixed stepping and restart" begin
    for iip in (false, true)
        f(u, p, t) = -u
        f!(du, u, p, t) = (du .= -u; nothing)
        prob = iip ? ODEProblem(f!, [1.0], (0.0, 2.0)) : ODEProblem(f, 1.0, (0.0, 2.0))
        integ = init(prob, FBDF(time_filter = true); adaptive = false, dt = 0.01)
        for _ in 1:8
            step!(integ)
        end
        @test integ.cache.iters_from_event == 8
        @test integ.cache.order > 1
        @test only(integ.u) ≈ exp(-integ.t) rtol = 1.0e-4
        derivative_discontinuity!(integ, true)
        step!(integ)
        @test integ.cache.order == 1
        @test integ.cache.filter_order == 0
        @test integ.cache.iters_from_event == 1
        for _ in 1:8
            step!(integ)
        end
        @test integ.cache.order > 1
        @test only(integ.u) ≈ exp(-integ.t) rtol = 2.0e-4
    end
end

@testset "FBDF filtering with rejected steps and callbacks" begin
    for iip in (false, true)
        f(u, p, t) = -u
        f!(du, u, p, t) = (du .= -u; nothing)
        prob = iip ? ODEProblem(f!, [1.0], (0.0, 2.0)) : ODEProblem(f, 1.0, (0.0, 2.0))
        function jump!(integ)
            if iip
                integ.u .*= 2
            else
                integ.u *= 2
            end
        end
        cb = DiscreteCallback((u, t, integ) -> t == 1, jump!; save_positions = (false, true))
        sol = solve(
            prob, FBDF(time_filter = true); callback = cb, tstops = [1.0],
            dt = 0.5, abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.stats.nreject > 0
        @test only(sol(0.5)) ≈ exp(-0.5) rtol = 1.0e-6
        @test only(sol(1.5)) ≈ 2exp(-1.5) rtol = 1.0e-6
        @test only(sol.u[end]) ≈ 2exp(-2) rtol = 1.0e-6
    end
end

@testset "FBDF filtered local convergence" begin
    for iip in (false, true), k in 1:4, spacing in (0.5, 1.0, 3.0, (0.5, 1.5, 0.75, 2.5, 1.25, 2.0))
        errors = map((0.2, 0.1)) do h
            integ = polynomial_history(k, h, spacing, iip; degree = k + 2)
            OrdinaryDiffEqBDF.perform_step!(integ, integ.cache)
            abs(only(integ.u) - h^(k + 2))
        end
        @test errors[1] / errors[2] ≈ 2^(k + 2) rtol = 1.0e-8
    end
end
