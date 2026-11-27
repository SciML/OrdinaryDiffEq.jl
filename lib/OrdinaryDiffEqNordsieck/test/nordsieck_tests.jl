using OrdinaryDiffEqNordsieck
using OrdinaryDiffEqNordsieck: JVODEController
using OrdinaryDiffEqCore: AbstractController, AbstractControllerCache, DummyController,
    DummyControllerCache, accept_step_controller, get_qmin, get_qmax, get_qmax_first_step,
    get_qsteady_min, get_qsteady_max
using SciMLBase: IntervalNonlinearProblem
import OrdinaryDiffEqCore: get_EEst, post_newton_controller!, set_EEst!, setup_controller_cache
using OrdinaryDiffEqTsit5: Tsit5
using DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_bigfloatlinear,
    prob_ode_bigfloat2Dlinear,
    prob_ode_linear, prob_ode_2Dlinear

probArr = [prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear]
testTol = 0.25
dts = 1 .// (2 .^ (10:-1:5))

@testset "Nordsieck Convergence Tests" begin
    for i in eachindex(probArr)
        sim = test_convergence(dts, probArr[i], AN5())
        @test sim.𝒪est[:final] ≈ 5 atol = testTol
        @test sim.𝒪est[:l2] ≈ 5 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 5 atol = testTol
    end
end

@testset "JVODE local error rejection recovery" begin
    function lorenz!(du, u, p, t)
        du[1] = 10 * (u[2] - u[1])
        du[2] = u[1] * (28 - u[3]) - u[2]
        du[3] = u[1] * u[2] - (8 / 3) * u[3]
        return nothing
    end
    lorenz(u, p, t) = (du = similar(u); lorenz!(du, u, p, t); du)
    for f in (lorenz!, lorenz)
        prob = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 2.0))
        ref = solve(prob, Tsit5(); reltol = 1.0e-12, abstol = 1.0e-12)
        for alg in (JVODE_Adams(), JVODE_BDF())
            sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6)
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 2.0
            @test sol.stats.nreject > 0
            @test norm(sol.u[end] - ref.u[end], Inf) < 1.0e-3
        end
    end
end

const JVODEFailureSnapshot = @NamedTuple{
    dt::Float64, force_stepfail::Bool, dts::Vector{Float64}, retry_dt::Float64,
}

struct JVODEHookController <: AbstractController
    snapshots::Vector{JVODEFailureSnapshot}
    quarter_step::Bool
end

struct JVODEHookControllerCache{C, S} <: AbstractControllerCache
    controller::C
    stock::S
    snapshots::Vector{JVODEFailureSnapshot}
    quarter_step::Bool
end

function setup_controller_cache(alg, cache, controller::JVODEHookController, ::Type{E}, disco_probs) where {E}
    stock = setup_controller_cache(alg, cache, JVODEController(; failfactor = 3), E, disco_probs)
    return JVODEHookControllerCache(stock.controller, stock, controller.snapshots, controller.quarter_step)
end

get_EEst(cache::JVODEHookControllerCache) = get_EEst(cache.stock)
set_EEst!(cache::JVODEHookControllerCache, val) = set_EEst!(cache.stock, val)

function post_newton_controller!(integrator, cache::JVODEHookControllerCache, alg)
    dt = integrator.dt
    force_stepfail = integrator.force_stepfail
    dts = copy(integrator.cache.dts)
    if cache.quarter_step
        integrator.dt /= 4
    else
        post_newton_controller!(integrator, cache.stock, alg)
    end
    push!(cache.snapshots, (; dt, force_stepfail, dts, retry_dt = integrator.dt))
    return nothing
end

@testset "JVODE corrector failure recovery" begin
    rhs(u, t) = -20 * (u - cos(t)) - sin(t)
    f(u, p, t) = [rhs(u[1], t)]
    function f!(du, u, p, t)
        du[1] = rhs(u[1], t)
        return nothing
    end
    opts = (; dt = 1.0, reltol = 1.0e-7, abstol = 1.0e-7)

    @testset "$(alg.algorithm), in-place = $(SciMLBase.isinplace(prob))" for
        alg in (JVODE_Adams(), JVODE_BDF()),
            prob in (ODEProblem(f!, [1.0], (0.0, 1.0)), ODEProblem(f, [1.0], (0.0, 1.0)))
        @testset "DummyController" begin
            @testset "Direct failure hook" begin
                alg_with_qmax = JVODE(alg.algorithm; qmax = 5.0)
                integrator = init(prob, alg_with_qmax; opts..., controller = DummyController())
                @test get_qmax(integrator) == 5
                @test integrator.dt == opts.dt
                post_newton_controller!(integrator, alg_with_qmax)
                @test integrator.dt == opts.dt / 2
            end

            @testset "Compact rational step-size options" begin
                bounds = map(
                    Rational{Int8},
                    (; qmin = 1 // 4, qmax = 5 // 1, qsteady_min = 9 // 8, qsteady_max = 7 // 4),
                )
                compact_alg = JVODE(alg.algorithm; bounds...)
                integrator = init(prob, compact_alg; opts..., controller = DummyController())
                getters = (get_qmin, get_qmax, get_qsteady_min, get_qsteady_max)
                resolved = map(getter -> getter(integrator), getters)
                @test resolved === values(bounds)
                @test integrator.dt == opts.dt
                post_newton_controller!(integrator, compact_alg)
                @test integrator.dt == opts.dt / 2

                sol = solve(prob, compact_alg; opts..., controller = DummyController())
                @test SciMLBase.successful_retcode(sol)
                @test sol.t[end] == 1.0
                @test only(sol.u[end]) ≈ cos(1) atol = 1.0e-6
            end

            sol = solve(prob, alg; opts..., controller = DummyController())
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 1.0
            @test only(sol.u[end]) ≈ cos(1) atol = 1.0e-6
        end

        @testset "Reuse a dedicated controller" begin
            bounds = (; qmin = 1 // 4, qmax = 5 // 1, qsteady_min = 9 // 8, qsteady_max = 7 // 4)
            controller = JVODEController(; bounds...)
            initial = init(prob, alg; opts..., controller)
            resolved = initial.controller_cache.controller
            event = ContinuousCallback(
                (u, t, integrator) -> t - 0.5, integrator -> nothing; maybe_discontinuity = true,
            )
            reused = init(prob, alg; opts..., controller = resolved, callback = event)
            getters = (get_qmin, get_qmax, get_qsteady_min, get_qsteady_max)
            @test map(getter -> getter(reused), getters) === values(bounds)
            @test isempty(resolved.basic.disco_probs)
            @test length(reused.controller_cache.controller.basic.disco_probs) == 1
            sol = solve!(reused)
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 1.0
            @test only(sol.u[end]) ≈ cos(1) atol = 1.0e-6
        end

        @testset "$name" for (name, quarter_step, factor) in (
                ("Custom post-Newton hook", true, 4),
                ("Configured failure factor", false, 3),
            )
            snapshots = JVODEFailureSnapshot[]
            controller = JVODEHookController(snapshots, quarter_step)
            integrator = init(prob, alg; opts..., controller)
            history = copy(integrator.cache.dts)
            # The first trial initializes the leading history entry before shifting it.
            history[1] = opts.dt
            step!(integrator)
            trials = first(snapshots, 2)
            @test [trial.dt for trial in trials] == [1.0, 1 / factor]
            for trial in trials
                @test trial.force_stepfail
                @test trial.dts == history
                @test trial.retry_dt == trial.dt / factor
            end
            sol = solve!(integrator)
            @test SciMLBase.successful_retcode(sol)
            @test sol.t[end] == 1.0
            @test only(sol.u[end]) ≈ cos(1) atol = 1.0e-6
        end
    end
end

@testset "JVODE DummyController retry with rational time" begin
    h = 1 // 2^26
    prob = ODEProblem((u, p, t) -> -u / h, 1.0, (0 // 1, h))
    bounds = map(Float16, (; qmin = 0.25, qmax = 5, qsteady_min = 0.25, qsteady_max = 1.5))
    alg = JVODE_BDF(; bounds...)
    tol = 0.1
    opts = (; dt = h, controller = DummyController(), abstol = tol, reltol = tol)

    integrator = init(prob, alg; opts...)
    post_newton_controller!(integrator, alg)
    @test integrator.dt == h / 2

    sol = solve(prob, alg; opts...)
    @test SciMLBase.successful_retcode(sol)
    @test sol.t[end] == h
    @test sol.u[end] ≈ exp(-1) atol = tol rtol = 0
end

@testset "JVODEController retry with rational time ($(alg.algorithm), flag = $flag)" for
    alg in (JVODE_Adams(), JVODE_BDF()), flag in (0.0, false)
    h = 1 // 2^26
    prob = ODEProblem((u, p, t) -> -u / h, 1.0, (0 // 1, h))
    tol = 0.1
    controller = JVODEController(; failfactor = Float16(2), discontinuity_detection = flag)
    opts = (; dt = h, controller, abstol = tol, reltol = tol)

    integrator = init(prob, alg; opts...)
    post_newton_controller!(integrator, alg)
    @test integrator.dt == h / 2

    sol = solve(prob, alg; opts...)
    @test SciMLBase.successful_retcode(sol)
    @test sol.t[2] == h / 2
    @test sol.t[end] == h
    @test sol.u[end] ≈ exp(-1) atol = tol rtol = 0
end

@testset "JVODEController rational retry ($(alg.algorithm), factor = $factor)" for
    alg in (JVODE_Adams(), JVODE_BDF()),
        (h, factor, retry, bounds) in (
            (1 // 10001, Float64(32768 // 10001), 1 // 32768, (; qsteady_min = 0.0, qsteady_max = Inf)),
            (9 // 5, 1.8, 1 // 1, NamedTuple()),
        )
    prob = ODEProblem((u, p, t) -> -u / h, 1.0, (0 // 1, h))
    tol = 0.1
    controller = JVODEController(; failfactor = factor, bounds...)
    opts = (; dt = h, controller, abstol = tol, reltol = tol)

    @testset "Direct failure hook" begin
        integrator = init(prob, alg; opts...)
        post_newton_controller!(integrator, alg)
        @test integrator.dt == retry
    end

    @testset "Solve" begin
        sol = solve(prob, alg; opts...)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t[2] == retry
        @test sol.t[end] == h
        @test sol.u[end] ≈ exp(-1) atol = tol rtol = 0
    end
end

@testset "JVODEController compact rational override ($(alg.algorithm), flag = $flag)" for
    alg in (JVODE_Adams(), JVODE_BDF()), flag in (0, false)
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    controller = JVODEController(; qmin = Rational{Int8}(1, 4), discontinuity_detection = flag)
    integrator = init(prob, alg; controller, reltol = 1.0e-8, abstol = 1.0e-8)
    @test get_qmin(integrator) === 1 // 4
    @test get_qmax_first_step(integrator) === 10000 // 1
    sol = solve!(integrator)
    @test SciMLBase.successful_retcode(sol)
    @test sol.t[end] == 1.0
    @test sol.u[end] ≈ exp(-1) atol = 1.0e-6
end

@testset "JVODE controller caches keep the requested error-estimate type" begin
    alg = JVODE_BDF()
    integrator = init(ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0)), alg)
    @testset "$(nameof(typeof(controller)))" for controller in (JVODEController(), DummyController())
        cache = setup_controller_cache(alg, integrator.cache, controller, AbstractFloat)
        set_EEst!(cache, nextfloat(BigFloat(1)))
        @test get_EEst(cache) isa BigFloat
        @test get_EEst(cache) > 1
        @test !accept_step_controller(integrator, cache, alg)

        cache = setup_controller_cache(alg, integrator.cache, controller, Real)
        set_EEst!(cache, 0.375)
        @test get_EEst(cache) === 0.375
    end
end

@testset "JVODE discontinuity options are independent of keyword order" begin
    prob = ODEProblem((u, p, t) -> -20 * (u - cos(t)) - sin(t), 1.0, (0.0, 1.0))
    alg = JVODE_BDF()
    overrides = (;
        qmin = 0.25, qmax = 5, qsteady_min = 0.25, qsteady_max = 1.5,
        qmax_first_step = 10000, gamma = 0.9, failfactor = 3,
    )
    @testset "$kind, disco_probs $order" for
        (kind, discontinuity_detection, disco_probs) in (
                ("integer options", 0, 0),
                ("typed options", false, IntervalNonlinearProblem[]),
            ),
            (order, kwargs) in (
                ("last", (; overrides..., discontinuity_detection, disco_probs)),
                ("first", (; disco_probs, overrides..., discontinuity_detection)),
            )
        events = Ref(0)
        callback = ContinuousCallback(
            (u, t, integrator) -> t - 0.5, integrator -> (events[] += 1);
            maybe_discontinuity = true,
        )
        integrator = init(
            prob, alg; controller = JVODEController(; kwargs...), callback,
            reltol = 1.0e-7, abstol = 1.0e-7,
        )
        basic = integrator.controller_cache.controller.basic
        @test map(key -> getproperty(basic, key), keys(overrides)) === values(map(Float64, overrides))
        @test basic.discontinuity_detection === false
        @test length(basic.disco_probs) == 1
        @test basic.disco_probs !== disco_probs
        sol = solve!(integrator)
        @test events[] == 1
        @test SciMLBase.successful_retcode(sol)
        @test sol.t[end] == 1.0
        @test sol.u[end] ≈ cos(1) atol = 1.0e-6
    end
end

struct JVODEViewController <: AbstractController end

checked_dummy_cache(cache::DummyControllerCache) = cache

function setup_controller_cache(alg, cache, ::JVODEViewController, ::Type{E}, problems) where {E}
    stock = setup_controller_cache(alg, cache, DummyController(), E, view(problems, eachindex(problems)))
    return checked_dummy_cache(stock)
end

struct JVODEWrappedCacheController <: AbstractController end

struct JVODEWrappedCache{C}
    inner::C
end

function setup_controller_cache(alg, cache::JVODEWrappedCache, controller::DummyController, ::Type{E}, problems) where {E}
    return setup_controller_cache(alg, cache.inner, controller, E, problems)
end

function setup_controller_cache(alg, cache, ::JVODEWrappedCacheController, ::Type{E}, problems) where {E}
    return setup_controller_cache(alg, JVODEWrappedCache(cache), DummyController(), E, problems)
end

@testset "JVODE DummyController cache delegation" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    @testset "$(alg.algorithm), $(nameof(typeof(controller)))" for
        alg in (JVODE_Adams(), JVODE_BDF()),
            controller in (JVODEViewController(), JVODEWrappedCacheController())
        events = Float64[]
        callback = ContinuousCallback(
            (u, t, integrator) -> t - 0.5, integrator -> push!(events, integrator.t);
            maybe_discontinuity = true,
        )
        integrator = init(
            prob, alg; controller, callback,
            reltol = 1.0e-8, abstol = 1.0e-8,
        )
        @test integrator.controller_cache isa DummyControllerCache
        @test integrator.controller_cache.cache === integrator.cache
        sol = solve!(integrator)
        @test length(events) == 1
        @test only(events) ≈ 0.5
        @test SciMLBase.successful_retcode(sol)
        @test sol.t[end] == 1.0
        @test sol.u[end] ≈ exp(-1) atol = 1.0e-6
    end
end

@testset "JVODE_BDF uses BDF formulas" begin
    oscillator!(du, u, p, t) = (du[1] = u[2]; du[2] = -u[1]; nothing)
    prob = ODEProblem(oscillator!, [1.0, 0.0], (0.0, 10.0))
    integ = init(prob, JVODE_BDF(); adaptive = false, dt = 1 / 64)
    for _ in 1:200
        step!(integ)
    end
    q = integ.cache.order
    @test q >= 2
    # Constant-step BDF coefficients in CVODE's normalization: the coefficients
    # of the polynomial ∏ⱼ₌₁^q (1 + x/j)
    expected = foldl((c, j) -> [c; 0] .+ [0; c] ./ j, 1:q; init = [1.0])
    @test integ.cache.l[1:(q + 1)] ≈ expected

    for prob in (prob_ode_linear, prob_ode_2Dlinear), tol in (1.0e-8, 1.0e-10)
        integ = init(prob, JVODE_BDF(); reltol = tol, abstol = tol)
        maxorder = 0
        while integ.t < prob.tspan[2]
            step!(integ)
            maxorder = max(maxorder, integ.cache.order)
        end
        @test maxorder == 5
    end

    sol_adams = solve(prob_ode_2Dlinear, JVODE_Adams(); reltol = 1.0e-6, abstol = 1.0e-6)
    sol_bdf = solve(prob_ode_2Dlinear, JVODE_BDF(); reltol = 1.0e-6, abstol = 1.0e-6)
    @test sol_adams.t != sol_bdf.t

    function rober!(du, u, p, t)
        y1, y2, y3 = u
        du[1] = -0.04y1 + 1.0e4 * y2 * y3
        du[2] = 0.04y1 - 3.0e7 * y2^2 - 1.0e4 * y2 * y3
        du[3] = 3.0e7 * y2^2
        return nothing
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0))
    # Rodas5P with reltol = 1e-13, abstol = 1e-14
    ref = [0.9664597373330047, 3.0746265785776e-5, 0.03350951640120991]
    sol_adams = solve(prob, JVODE_Adams(); reltol = 1.0e-6, abstol = 1.0e-6)
    sol_bdf = solve(prob, JVODE_BDF(); reltol = 1.0e-6, abstol = 1.0e-6)
    @test SciMLBase.successful_retcode(sol_bdf)
    @test sol_bdf.stats.naccept < sol_adams.stats.naccept
    @test norm(sol_bdf.u[end] - ref, Inf) < 1.0e-5
end

probArr = [
    prob_ode_linear,
    prob_ode_2Dlinear,
]
@testset "Nordsieck Adaptivity Tests: AN5" begin
    for i in eachindex(probArr)
        prob = probArr[i]
        sol = solve(prob, AN5(), reltol = 1.0e-6)
        @test length(sol.t) < 11
        @test SciMLBase.successful_retcode(sol)
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test exact ≈ sol.u[end] atol = 1.0e-5
    end
end

@testset "Nordsieck Adaptivity Tests: JVODE" begin
    for i in eachindex(probArr),
            sol in [JVODE_Adams(), JVODE_BDF()]

        prob = probArr[i]
        sol = solve(prob, sol, reltol = 1.0e-4, abstol = 1.0e-7)
        @test length(sol.t) < 22
        @test SciMLBase.successful_retcode(sol)
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test norm(exact - sol.u[end], Inf) < 3.0e-3
    end
end

const A_resize = [-1.0 0.35 0.1; 0.2 -1.4 0.25; 0.15 0.3 -1.2]
const U0_resize = [1.0, -0.5]
const UNEW_resize = 0.3
const TG_resize = 0.5
const TS_resize = 0.65
const TF_resize = 0.8

f_resize!(du, u, p, t) = (mul!(du, view(A_resize, 1:length(u), 1:length(u)), u); nothing)

function exact_resize(t)
    u = exp(A_resize[1:2, 1:2] * min(t, TG_resize)) * U0_resize
    t <= TG_resize && return u
    u3 = exp(A_resize * (min(t, TS_resize) - TG_resize)) * vcat(u, UNEW_resize)
    t <= TS_resize && return u3
    return exp(A_resize[1:2, 1:2] * (t - TS_resize)) * u3[1:2]
end

resize_cases = (
    (AN5(), 1.0e-10, 1.0e-12, 1.0e-7),
    (JVODE_Adams(), 1.0e-8, 1.0e-10, 1.0e-6),
    (JVODE_BDF(), 1.0e-8, 1.0e-10, 1.0e-6),
)
@testset "Nordsieck state resizing ($alg, u_modified=$modified)" for
    (alg, reltol, abstol, check_rtol) in resize_cases, modified in (true, false)
    grew = Ref(false)
    shrank = Ref(false)
    restart_orders = Int[]
    function affect!(integrator)
        n = grew[] ? 2 : 3
        resize!(integrator, n)
        n == 3 && (integrator.u[3] = UNEW_resize)
        u_modified!(integrator, modified)
        push!(restart_orders, integrator.cache.order)
        grew[] ? (shrank[] = true) : (grew[] = true)
        return nothing
    end
    cb = DiscreteCallback(
        (u, t, integrator) -> (!grew[] && t >= TG_resize) || (!shrank[] && t >= TS_resize),
        affect!;
        save_positions = (false, false)
    )
    sol = solve(
        ODEProblem(f_resize!, copy(U0_resize), (0.0, TF_resize)), alg;
        callback = cb, tstops = [TG_resize, TS_resize], reltol = reltol, abstol = abstol
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol(0.58) ≈ exact_resize(0.58) rtol = check_rtol atol = 1.0e-10
    @test sol.u[end] ≈ exact_resize(TF_resize) rtol = check_rtol atol = 1.0e-10
    @test restart_orders == [1, 1]
end

@testset "resize before first step initializes restart state" begin
    for alg in (AN5(), JVODE_Adams(), JVODE_BDF())
        it = init(ODEProblem(f_resize!, [1.0, -0.5], (0.0, 0.1)), alg; adaptive = false, dt = 1.0e-4)
        resize!(it, 3)
        it.u[3] = 0.3
        it.uprev[3] = NaN
        u_modified!(it, true)
        step!(it)
        @test all(isfinite, it.u)
    end
end
