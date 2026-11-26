using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, Test

state_value(u::Number) = u
state_value(u) = only(u)

function startup_problem(rhs, u0, tspan, iip)
    return if iip
        ODEProblem((du, u, p, t) -> (du[1] = rhs(u[1], t)), [u0], tspan)
    else
        ODEProblem((u, p, t) -> rhs(u, t), u0, tspan)
    end
end

function jump_callback(iip, jump)
    affect! = iip ? (integ -> (integ.u[1] += jump)) : (integ -> (integ.u += jump))
    return DiscreteCallback((u, t, integ) -> t == 1, affect!)
end

@testset "FBDF restart at a time stop" begin
    for iip in (false, true), filter in (false, true)
        @testset "iip=$iip filter=$filter" begin
            alg = FBDF(time_filter = filter)
            callback = jump_callback(iip, 0.25)
            rhs(u, t) = 2t
            abstol, reltol = 1.0e-10, 1.0e-8

            prob = startup_problem(rhs, 1.0, (0.0, 1.0003), iip)
            integ = init(prob, alg; callback, tstops = [1.0], dt = 0.01, abstol, reltol)
            while integ.t < 1
                step!(integ)
            end
            before = state_value(integ.u)
            step!(integ)
            increment = state_value(integ.u) - before
            @test integ.t > 1
            @test abs(increment - (integ.t^2 - 1)) <= abstol + reltol * abs(before)

            prob = startup_problem(rhs, 1.0, (0.0, 1.02), iip)
            integ = init(prob, alg; callback, tstops = [1.0], dt = 0.01, adaptive = false, abstol, reltol)
            while integ.t < 1
                step!(integ)
            end
            step!(integ)
            @test integ.t == 1.0 + 0.01
        end
    end
end

@testset "FBDF SDIRK restart retries" begin
    for iip in (false, true)
        @testset "iip=$iip" begin
            after_jump = Ref(false)
            rejected = Ref(false)
            trials = Real[]
            rhs(u, t) = begin
                after_jump[] || return zero(u)
                push!(trials, u)
                -u
            end
            affect! = function (integ)
                if iip
                    integ.u[1] = 7.0e-5
                else
                    integ.u = 7.0e-5
                end
                after_jump[] = true
                return
            end
            callback = DiscreteCallback((u, t, integ) -> t == 1, affect!)
            isoutofdomain = function (u, p, t)
                if t > 1 && !rejected[]
                    rejected[] = true
                    return true
                end
                return false
            end
            prob = startup_problem(rhs, 1.0, (0.0, 1.1), iip)
            abstol, reltol = 1.0e-11, 1.0e-6
            # The stop just before the jump gives the previous interpolant a
            # very short interval compared with the forced restart attempt.
            tstops = [1.0 - 1.0e-12, 1.0]
            integ = init(prob, FBDF(rk_start = true); callback, tstops, isoutofdomain, abstol, reltol)
            while integ.t < 1
                step!(integ)
            end
            set_proposed_dt!(integ, 0.05)
            step!(integ)
            @test rejected[]
            @test !isempty(trials)
            @test 6.0e-5 <= minimum(trials) <= maximum(trials) <= 8.0e-5
            @test integ.t > 1
            @test state_value(integ.u) ≈ 7.0e-5 * exp(1 - integ.t) atol = abstol rtol = reltol
        end
    end
end

@testset "FBDF SDIRK startup" begin
    for iip in (false, true)
        @testset "iip=$iip" begin
            prob = startup_problem((u, t) -> -u, 1.0, (0.0, 0.1), iip)
            integ = init(prob, FBDF(rk_start = true); abstol = 1.0e-9, reltol = 1.0e-7)
            step!(integ)
            @test integ.t > 0
            @test state_value(integ.u) ≈ exp(-integ.t) rtol = 1.0e-7
            step!(integ)
            @test integ.cache.prev_order == 2
            @test integ.cache.order == 2

            # For u' = -u and u(0) = 1, the scaled norms of f₀ and
            # (f₁ - f₀) / h₀ in the Hairer–Wanner proposal are both 1 / scale.
            scale = 1.0e-8 + 1.0e-6
            h = min(100 * 0.01, sqrt(0.01 * scale), 0.1)
            integ = init(prob, FBDF(); abstol = 1.0e-8, reltol = 1.0e-6)
            @test integ.dt ≈ h rtol = 100eps()

            errors = map((0.02, 0.01, 0.005)) do dt
                integ = init(prob, FBDF(rk_start = true); dt, abstol = 1.0e-3, reltol = 1.0e-3)
                step!(integ)
                @test integ.t == dt
                return (
                    abs(state_value(integ.u) - exp(-dt)),
                    abs(state_value(integ.sol(dt / 2)) - exp(-dt / 2)),
                )
            end
            for j in 1:2, k in 1:2
                @test errors[k][j] / errors[k + 1][j] > 7
            end

            quadratic = startup_problem((u, t) -> 2t, 1.0, (0.0, 0.1), iip)
            integ = init(quadratic, FBDF(rk_start = true); dt = 0.01, abstol = 1.0e-4, reltol = 1.0e-4)
            step!(integ)
            @test integ.t == 0.01
            @test state_value(integ.sol(0.005)) ≈ 1 + 0.005^2 atol = 8eps() rtol = 0
        end
    end
    @test haskey(Base.Docs.meta(OrdinaryDiffEqCore), Base.Docs.Binding(OrdinaryDiffEqCore, :initial_timestep_order))
end

@testset "FBDF SDIRK startup resolves stiff decay" begin
    λ = -1.0e12
    dt = 0.01
    abstol, reltol = 1.0e-8, 1.0e-6
    tspan = (0.0, dt)
    alg = FBDF(rk_start = true)

    for iip in (false, true)
        @testset "iip=$iip" begin
            prob = startup_problem((u, t) -> λ * u, 1.0, tspan, iip)
            integ = init(prob, alg; dt, abstol, reltol)
            step!(integ)
            tmid = integ.t / 2
            exact = exp(λ * tmid)
            @test abs(state_value(integ.sol(tmid)) - exact) <= abstol + reltol * abs(exact)

            fired = Ref(false)
            function terminate_on_crossing!(integ)
                fired[] = true
                return terminate!(integ)
            end
            callback = ContinuousCallback(
                (u, t, integ) -> state_value(u) + 0.5, nothing, terminate_on_crossing!
            )
            sol = solve(prob, alg; callback, dt, abstol, reltol)
            @test !fired[]
            @test sol.retcode == ReturnCode.Success
            @test sol.t[end] == last(tspan)
        end
    end
end
