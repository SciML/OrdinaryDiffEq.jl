using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqRosenbrock, Test

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/4801
# u1(t) = 1/(1 - t) exactly; u2 is a clock (u2 == t for any consistent method).
f_blowup!(du, u, p, t) = (du[1] = u[1]^2; du[2] = 1.0; nothing)
f_blowup(u, p, t) = [u[1]^2, 1.0]

function f_rober!(du, u, p, t)
    y1, y2, y3 = u
    du[1] = -0.04y1 + 1.0e4 * y2 * y3
    du[2] = 0.04y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2^2
    du[3] = 3.0e7 * y2^2
    return nothing
end
f_rober(u, p, t) = (du = similar(u); f_rober!(du, u, p, t); du)

blowup_probs = (
    ODEProblem(f_blowup!, [1.0, 0.0], (0.0, 0.9)),
    ODEProblem(f_blowup, [1.0, 0.0], (0.0, 0.9)),
)
rober_probs = (
    ODEProblem(f_rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5)),
    ODEProblem(f_rober, [1.0, 0.0, 0.0], (0.0, 1.0e5)),
)
rober_ref = solve(rober_probs[1], Rodas5P(); reltol = 1.0e-12, abstol = 1.0e-14).u[end]

pid_controllers(alg) = (
    OrdinaryDiffEqCore.PIDController(0.7, -0.4),
    OrdinaryDiffEqCore.PIDController(alg),
)

# A run is stalled if it needs far more steps than the default controller (master took
# ~1e6 accepts at dt ~ 1e-8 where the default controller takes ~1e3).
const MAXITERS = 100_000

@testset "PIDController with BDF-family methods" begin
    @testset "PIDController acceptance agrees with the step commit" begin
        for alg in (QNDF1(), QNDF2(), QNDF(), FBDF(), ABDF2()), prob in blowup_probs,
                controller in pid_controllers(alg)
            integ = init(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6, controller)
            fixed_order = alg isa Union{QNDF1, QNDF2, ABDF2}
            mismatched = 0
            while integ.t < prob.tspan[2] && integ.iter < MAXITERS
                step!(integ)
                integ.accept_step || continue
                # Fixed-order methods follow PID's own acceptance rule and must record
                # every accepted step; variable-order methods accept on the error estimate.
                if fixed_order
                    integ.cache.dtₙ₋₁ == integ.dt || (mismatched += 1)
                else
                    OrdinaryDiffEqCore.get_EEst(integ) > 1 && (mismatched += 1)
                end
            end
            @test mismatched == 0
            @test integ.t == prob.tspan[2]
        end
    end

    @testset "$(nameof(typeof(alg))) + PIDController on the blowup problem" for alg in (
            QNDF1(), QNDF2(), QNDF(), FBDF(), ABDF2(),
        )
        for prob in blowup_probs
            default = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6)
            for controller in pid_controllers(alg)
                sol = solve(
                    prob, alg; reltol = 1.0e-6, abstol = 1.0e-6, controller, maxiters = MAXITERS
                )
                @test sol.retcode == ReturnCode.Success
                @test sol.t[end] == 0.9
                # The clock is exact up to the Newton tolerance, so a lost step shows here.
                @test sol.u[end][2] ≈ 0.9 atol = 1.0e-6
                # Exact value 1/(1 - 0.9) = 10; the default controller is within 1.1% here.
                @test sol.u[end][1] ≈ 10 rtol = 0.02
                @test sol.stats.naccept <= 2 * default.stats.naccept
            end
        end
    end

    @testset "$(nameof(typeof(alg))) + PIDController on Robertson" for alg in (
            QNDF1(), QNDF2(), QNDF(), FBDF(),
        )
        for prob in rober_probs
            default = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-10)
            for controller in pid_controllers(alg)
                sol = solve(
                    prob, alg; reltol = 1.0e-6, abstol = 1.0e-10, controller, maxiters = MAXITERS
                )
                @test sol.retcode == ReturnCode.Success
                @test sol.t[end] == 1.0e5
                # The default controller reaches y1, y3 within 9e-4 relative of the Rodas5P
                # reference with QNDF1, the least accurate method here.
                @test sol.u[end][1] ≈ rober_ref[1] rtol = 5.0e-3
                @test sol.u[end][3] ≈ rober_ref[3] rtol = 5.0e-3
                @test sum(sol.u[end]) ≈ 1 atol = 1.0e-8
                @test sol.stats.naccept <= 2 * default.stats.naccept
            end
        end
    end
end

struct VetoController <: OrdinaryDiffEqCore.AbstractController
    calls::Base.RefValue{Int}
    vetoes::Base.RefValue{Int}
    rejects_at_veto::Base.RefValue{Int}
    transition::Vector{Tuple{Float64, Int}}
end

mutable struct VetoCache{C, E} <: OrdinaryDiffEqCore.AbstractControllerCache
    controller::C
    EEst::E
    calls::Base.RefValue{Int}
    vetoes::Base.RefValue{Int}
    rejects_at_veto::Base.RefValue{Int}
    transition::Vector{Tuple{Float64, Int}}
end

function OrdinaryDiffEqCore.setup_controller_cache(alg, _, c::VetoController, ::Type{E}, disco_probs) where {E}
    basic = OrdinaryDiffEqCore.resolve_basic((;), alg, Float64; disco_probs)
    return VetoCache((; basic), oneunit(E), c.calls, c.vetoes, c.rejects_at_veto, c.transition)
end

function OrdinaryDiffEqCore.accept_step_controller(integrator, cache::VetoCache, alg)
    cache.calls[] += 1
    if length(cache.transition) == 1
        push!(cache.transition, (integrator.t, integrator.stats.nreject))
    end
    if OrdinaryDiffEqCore.get_EEst(integrator) <= 1 && cache.vetoes[] == 0
        push!(cache.transition, (integrator.t, integrator.stats.nreject))
        cache.vetoes[] += 1
        cache.rejects_at_veto[] = integrator.stats.nreject
        return false
    end
    return OrdinaryDiffEqCore.get_EEst(integrator) <= 1
end

@testset "BDF controller acceptance regressions" begin
    calls, vetoes, rejects_at_veto = Ref(0), Ref(0), Ref(0)
    transition = Tuple{Float64, Int}[]
    sol = solve(blowup_probs[2], FBDF(); controller = VetoController(calls, vetoes, rejects_at_veto, transition), abstol = 1.0e-6, reltol = 1.0e-6)
    @test calls[] > 0
    @test vetoes[] == 1
    @test sol.stats.nreject > rejects_at_veto[]
    @test length(transition) == 2
    @test transition[2][1] == transition[1][1]
    @test transition[2][2] == transition[1][2] + 1
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end][2] ≈ 0.9 atol = 1.0e-6

    controller = OrdinaryDiffEqCore.PIDController(1 // 6, 1 // 6; accept_safety = 1.0, retry_cap = 0.81)
    integ = init(blowup_probs[2], QNDF1(); controller, abstol = 1.0e-6, reltol = 1.0e-6)
    integ.controller_cache.dt_factor = 1.0
    old_dt = integ.dt
    OrdinaryDiffEqCore.step_reject_controller!(integ, integ.controller_cache, integ.alg)
    @test 0 < integ.dt < old_dt
    sol = solve(blowup_probs[2], QNDF1(); controller, abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 20_000)
    @test sol.stats.nreject > 0
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end][2] ≈ 0.9 atol = 1.0e-6

    f(u, p, t) = t^2
    f!(du, u, p, t) = (du[1] = t^2; nothing)
    for prob in (ODEProblem(f, 0.0, (0.0, 0.5)), ODEProblem(f!, [0.0], (0.0, 0.5)))
        integ = init(prob, ABDF2(); controller = OrdinaryDiffEqCore.PIDController(0.7, -0.4), dt = 0.125, dtmax = 0.125, abstol = 3 / 8192, reltol = 0.0, tstops = [0.1875])
        step!(integ)
        step!(integ)
        @test integ.t == 0.1875
        @test OrdinaryDiffEqCore.get_EEst(integ) == 1.0
        @test integ.cache.dtₙ₋₁ == 0.0625
        @test integ.cache.fsalfirstprev == (integ.u isa Number ? 0.125^2 : [0.125^2])
    end

    # The same tie rejected by `isoutofdomain` must leave the ABDF2 history at the
    # Euler step, and the retry must use that history.
    for prob in (ODEProblem(f, 0.0, (0.0, 0.5)), ODEProblem(f!, [0.0], (0.0, 0.5)))
        integ_ref = Ref{Any}()
        at_tie = Ref{Any}(nothing)
        isout = function (u, p, t)
            if t == 0.1875 && at_tie[] === nothing
                c = integ_ref[].cache
                at_tie[] = (OrdinaryDiffEqCore.get_EEst(integ_ref[]), c.dtₙ₋₁, first(c.fsalfirstprev))
            end
            return t == 0.1875 && first(u) > 0.0036
        end
        integ = init(prob, ABDF2(); dt = 0.125, dtmax = 0.125, abstol = 3 / 8192, reltol = 0.0, tstops = [0.1875], isoutofdomain = isout)
        integ_ref[] = integ
        step!(integ)
        u₁ = first(integ.u)
        step!(integ)
        @test at_tie[] == (1.0, 0.125, 0.0)
        @test integ.stats.nreject == 1
        h = integ.dt
        @test 0 < h < 0.0625
        @test integ.t ≈ 0.125 + h
        ρ = h / 0.125
        abdf2 = (1 + ρ^2 / 3) * u₁ + h * (2 / 3 * integ.t^2 - (ρ - 1) / 3 * 0.125^2)
        @test first(integ.u) ≈ abdf2 rtol = 1.0e-12
        @test integ.cache.dtₙ₋₁ == h
    end

    # `PIDController` accepts on its predicted factor: a limiter that returns 0.5 on its
    # first call rejects a step whose error estimate passes, and the history must follow.
    for alg in (QNDF1(), QNDF2()), iip in (false, true)
        calls = Ref(0)
        limiter(x) = (calls[] += 1; calls[] == 1 ? 0.5 : 1.0)
        prob = iip ? ODEProblem((du, u, p, t) -> (du[1] = 1.0; nothing), [0.0], (0.0, 0.5)) :
            ODEProblem((u, p, t) -> 1.0, 0.0, (0.0, 0.5))
        controller = OrdinaryDiffEqCore.PIDController(0.7, -0.4; limiter, accept_safety = 0.81)
        sol = solve(prob, alg; dt = 0.01, abstol = 1.0e-4, reltol = 1.0e-4, controller)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject == 1
        @test sol.t[2] == 0.0025
        @test first(sol.u[end]) ≈ 0.5 atol = 1.0e-12
    end
end

@testset "Six-member composite commits the active BDF history" begin
    for iip in (false, true)
        prob = iip ? ODEProblem((du, u, p, t) -> (du[1] = t^2; nothing), [0.0], (0.0, 1.0)) :
            ODEProblem((u, p, t) -> t^2, 0.0, (0.0, 1.0))
        ref = solve(prob, QNDF2(); dt = 1 / 8, adaptive = false)
        alg = OrdinaryDiffEqCore.CompositeAlgorithm(ntuple(_ -> QNDF2(), 6), _ -> 6)
        integ = init(prob, alg; dt = 1 / 8, adaptive = false)
        @test integ.cache isa OrdinaryDiffEqCore.DefaultCache
        recorded = Float64[]
        while integ.t < prob.tspan[2]
            step!(integ)
            @test integ.cache.current == 6
            push!(recorded, integ.cache.cache6.dtₙ₋₁)
        end
        @test recorded == fill(1 / 8, 8)
        @test integ.sol.t == ref.t
        @test first(integ.sol.u[end]) == first(ref.u[end])
    end
end

# A limiter that returns 7/8 once rejects the first step (accept_safety = 15/16). By
# default the retry uses that predicted factor; `retry_cap` bounds it from above.
@testset "PIDController retry length with a rejecting limiter" begin
    dt = 1 / 64
    for alg in (QNDF1(), QNDF2(), ABDF2()), iip in (false, true),
            (retry_cap, factor) in ((nothing, 7 / 8), (0.81, 0.81))
        calls = Ref(0)
        limiter(x) = (calls[] += 1; calls[] == 1 ? 7 / 8 : 1.0)
        prob = iip ? ODEProblem((du, u, p, t) -> (du[1] = 0.0; nothing), [0.0], (0.0, 0.5)) :
            ODEProblem((u, p, t) -> 0.0, 0.0, (0.0, 0.5))
        kwargs = retry_cap === nothing ? (;) : (; retry_cap)
        controller = OrdinaryDiffEqCore.PIDController(0.7, -0.4; limiter, accept_safety = 15 / 16, kwargs...)
        sol = solve(prob, alg; dt, abstol = 1.0e-4, reltol = 1.0e-4, controller)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject == 1
        @test sol.t[2] == dt * factor
        @test all(u -> iszero(first(u)), sol.u)
    end
    for retry_cap in (0, 1, -0.5, 1.5)
        @test_throws ArgumentError OrdinaryDiffEqCore.PIDController(0.7, -0.4; retry_cap)
    end
end

@testset "PIDController retry_cap precision and precedence" begin
    PID = OrdinaryDiffEqCore.PIDController
    oneshot(value) = (calls = Ref(0); x -> (calls[] += 1; calls[] == 1 ? value : 1.0))
    zero_rhs(iip) = iip ? ODEProblem((du, u, p, t) -> (du[1] = 0.0; nothing), [0.0], (0.0, 0.5)) :
        ODEProblem((u, p, t) -> 0.0, 0.0, (0.0, 0.5))
    dt = 1 / 8

    for cap in (prevfloat(1.0), nextfloat(0.0))
        @test_throws ArgumentError PID(Float32, QNDF1(); retry_cap = cap)
        @test_throws ArgumentError PID(0.7f0, -0.4f0; retry_cap = cap)
        controller = PID(0.7, -0.4; retry_cap = cap, qmin = 0.125f0)
        @test controller.retry_cap === cap
        @test_throws ArgumentError init(zero_rhs(false), QNDF1(); dt, controller)
    end
    @test PID(Float32, QNDF1(); retry_cap = 0.5).retry_cap === 0.5f0

    for alg in (QNDF1(), ABDF2()), iip in (false, true),
            (qmin, cap) in ((0.125f0, 0.5), (0.2, 1 / 8)), value in (NaN, -1.0)
        controller = PID(0.7, -0.4; limiter = oneshot(value), retry_cap = cap, qmin)
        integ = init(zero_rhs(iip), alg; dt, controller)
        @test integ.controller_cache.controller.retry_cap === typeof(qmin)(cap)
        sol = solve!(integ)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject == 1
        @test sol.t[2] == dt * cap
    end

    # Domain rejection retries with `qmin` and the solver clamps the retry to `dtmin`,
    # so both can exceed `dt * retry_cap`.
    for iip in (false, true)
        checks = Ref(0)
        isoutofdomain = (u, p, t) -> (checks[] += 1; checks[] == 1)
        controller = PID(0.7, -0.4; retry_cap = 1 / 8, qmin = 1 / 2)
        sol = solve(zero_rhs(iip), QNDF1(); dt, isoutofdomain, controller)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject == 1
        @test sol.t[2] == dt / 2

        controller = PID(0.7, -0.4; limiter = oneshot(NaN), retry_cap = 1 / 8)
        sol = solve(zero_rhs(iip), QNDF1(); dt, dtmin = 1 / 16, force_dtmin = true, controller)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject == 1
        @test sol.t[2] == 1 / 16
    end
end
