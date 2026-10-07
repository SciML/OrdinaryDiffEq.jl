using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, ForwardDiff, Test, LinearAlgebra
using SciMLBase: CheckInit, NoInit
using OrdinaryDiffEqCore: DEVerbosity
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, NLNewton
using RecursiveArrayTools: ArrayPartition

foop = (u, p, t) -> u * p
proboop = ODEProblem(foop, ones(2), (0.0, 1000.0), 1.0)

fiip = (du, u, p, t) -> du .= u .* p
probiip = ODEProblem(fiip, ones(2), (0.0, 1000.0), 1.0)

@testset "FBDF reinit" begin
    for prob in [proboop, probiip]
        integ = init(prob, FBDF(), verbose = DEVerbosity(SciMLLogging.None())) #suppress warning to clean up CI
        solve!(integ)
        @test integ.sol.retcode != ReturnCode.Success
        @test integ.sol.t[end] >= 700
        reinit!(integ, prob.u0)
        solve!(integ)
        @test integ.sol.retcode != ReturnCode.Success
        @test integ.sol.t[end] >= 700
    end
end

function ad_helper(alg, prob)
    return function costoop(p)
        _oprob = remake(prob; p)
        sol = solve(_oprob, alg, saveat = 1:10)
        return sum(sum, sol.u)
    end
end

@testset "parameter autodiff" begin
    for prob in [proboop, probiip]
        for alg in [FBDF(), QNDF()]
            ForwardDiff.derivative(ad_helper(alg, prob), 1.0)
        end
    end
end

@testset "FBDF with non-default max_order" begin
    # Test that FBDF works with max_order < 5 (regression test for hardcoded Val(5))
    # MO=1 (backward Euler) is more conservative with the CVODE step size formula,
    # so it reaches fewer time steps on exponential growth problems.
    for MO in 1:5
        for prob in [proboop, probiip]
            sol = solve(prob, FBDF(max_order = Val{MO}()), verbose = DEVerbosity(SciMLLogging.None()))
            @test sol.t[end] >= (MO == 1 ? 250 : 700)
        end
    end
end

@testset "DFBDF with non-default max_order" begin
    # Test that DFBDF works with max_order < 5 (regression test for hardcoded Val(5))
    # The bug caused BoundsError when max_order != 5 due to data/type mismatch
    function dfbdf_dae_f!(resid, du, u, p, t)
        resid[1] = -0.5 * u[1] + u[2] - du[1]
        resid[2] = u[1] - u[2] - du[2]
    end
    dae_prob_mo = DAEProblem(
        dfbdf_dae_f!, zeros(2), [1.0, 1.0], (0.0, 1.0),
        differential_vars = [true, true]
    )
    for MO in 1:5
        sol = solve(
            dae_prob_mo, DFBDF(max_order = Val{MO}()), initializealg = BrownFullBasicInit(),
            abstol = 1.0e-8, reltol = 1.0e-8, verbose = DEVerbosity(SciMLLogging.None())
        )
        @test sol.t[end] > 0
    end
end

# [sources] in Project.toml requires Julia ≥ 1.11 so the local OrdinaryDiffEqCore
# fix is only picked up in CI on 1.11+.
if VERSION >= v"1.11"
    @testset "get_du during init callback (issue #3117)" begin
        # get_du must not crash when called before the first step, e.g. from a
        # callback that fires during init. Before the fix, k was empty at that
        # point and the stiff interpolation threw a BoundsError.
        f_ode!(du, u, p, t) = (du[1] = -u[1])
        prob_ode = ODEProblem(f_ode!, [1.0], (0.0, 1.0))

        du_at_init = Ref{Vector{Float64}}()
        function init_cb(c, u, t, integrator)
            du_at_init[] = get_du(integrator)
        end
        cb = DiscreteCallback(
            (u, t, integrator) -> false, identity;
            initialize = init_cb
        )

        for alg in (QNDF(), FBDF())
            du_at_init[] = Float64[]
            integrator = init(prob_ode, alg; callback = cb, save_everystep = false)
            @test du_at_init[][1] ≈ -1.0
            solve!(integrator)
        end

        # Also test get_du! (in-place variant)
        du_buf = [0.0]
        function init_cb_ip(c, u, t, integrator)
            get_du!(du_buf, integrator)
        end
        cb_ip = DiscreteCallback(
            (u, t, integrator) -> false, identity;
            initialize = init_cb_ip
        )

        for alg in (QNDF(), FBDF())
            du_buf[1] = 0.0
            integrator = init(prob_ode, alg; callback = cb_ip, save_everystep = false)
            @test du_buf[1] ≈ -1.0
            solve!(integrator)
        end

        # DFBDF (DAE): get_du before first step should throw a clear error
        # because integrator.du is not initialized until the solver steps.
        function dae_f!(resid, du, u, p, t)
            resid[1] = du[1] + u[1]
        end
        dae_prob = DAEProblem(
            dae_f!, [-1.0], [1.0], (0.0, 1.0);
            differential_vars = [true]
        )

        dae_errored = Ref(false)
        function init_cb_dae(c, u, t, integrator)
            try
                get_du(integrator)
            catch e
                dae_errored[] = isa(e, ErrorException) &&
                    contains(e.msg, "DAE problems")
            end
        end
        cb_dae = DiscreteCallback(
            (u, t, integrator) -> false, identity;
            initialize = init_cb_dae
        )

        integrator = init(dae_prob, DFBDF(); callback = cb_dae, save_everystep = false)
        @test dae_errored[]
        solve!(integrator)
    end
end

if VERSION >= v"1.12"
    @testset "FBDF in-place perform_step! non-allocating" begin
        integrator = init(
            probiip, FBDF(), abstol = 1.0e-8, reltol = 1.0e-8,
            save_everystep = false
        )
        # Warm up to reach higher orders and compile all code paths
        for _ in 1:10
            step!(integrator)
        end
        allocs = @allocated step!(integrator)
        @test allocs == 0
    end

    @testset "DFBDF in-place perform_step! non-allocating" begin
        function dae_f!(resid, du, u, p, t)
            resid[1] = -0.5 * u[1] + u[2] - du[1]
            resid[2] = u[1] - u[2] - du[2]
        end
        dae_prob = DAEProblem(
            dae_f!, zeros(2), [1.0, 1.0], (0.0, 1.0),
            differential_vars = [true, false]
        )
        integrator = init(
            dae_prob, DFBDF(), abstol = 1.0e-8, reltol = 1.0e-8,
            save_everystep = false, initializealg = BrownFullBasicInit()
        )
        for _ in 1:10
            step!(integrator)
        end
        allocs = @allocated step!(integrator)
        @test allocs == 0
    end
end

# Regression test for issue #3645: Newton failure with QNDF used to recurse
# through `post_newton_controller!(integrator, alg::QNDF)` →
# generic 2-arg in core → `BDFControllerCache` 3-arg → 2-arg again …,
# producing a `StackOverflowError` on the first failed Newton step.
@testset "QNDF Newton failure does not StackOverflow (#3645)" begin
    f_qndf!(du, u, p, t) = (du[1] = -1.0e6 * (u[1] - cos(t)); nothing)
    prob_qndf = ODEProblem(f_qndf!, [0.0], (0.0, 1.0))
    # Cripple the Newton solver so it can never converge.
    alg = QNDF(; nlsolve = NLNewton(; max_iter = 1, κ = 1.0e-30))
    sol = solve(
        prob_qndf, alg; dt = 0.5, reltol = 1.0e-12, abstol = 1.0e-12,
        verbose = DEVerbosity(SciMLLogging.None())
    )
    # With the fix, repeated Newton failures shrink dt until dtmin and
    # the solver gives up cleanly (Unstable) instead of overflowing the stack.
    @test sol.retcode != ReturnCode.Success
    @test sol.retcode != ReturnCode.Default
end

# Regression test for issue #3962: the out-of-place Newton step used
# `_reshape(W \ _vec(ztmp), axes(ztmp))`, which collapses the ArrayPartition state
# behind SecondOrderODEProblem/DynamicalODEFunction down to a bare Vector. Storing
# `z .- dz` back into the ArrayPartition-typed nlsolver fields then failed with
# `Cannot convert Vector to ArrayPartition`.
@testset "OOP SecondOrderODEProblem with FBDF preserves ArrayPartition (#3962)" begin
    # u'' = -u  ⇒  position = [cos t, sin t], velocity = [-sin t, cos t]
    ho_iip(ddu, du, u, p, t) = (@. ddu = -u)
    ho_oop(du, u, p, t) = -u
    u0 = [1.0, 0.0]
    du0 = [0.0, 1.0]
    tspan = (0.0, 1.0)
    prob_iip = SecondOrderODEProblem(ho_iip, du0, u0, tspan)
    prob_oop = SecondOrderODEProblem(ho_oop, du0, u0, tspan)

    ref = solve(prob_iip, FBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
    sol = solve(prob_oop, FBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end] isa ArrayPartition
    @test norm(sol.u[end] - ref.u[end]) < 1.0e-6
end

@testset "backward-in-time step rejection shrinks |dt| (#4504)" begin
    # Signed step comparisons must compare magnitudes: for tdir < 0 both h and
    # hₖ₋₁ are negative, so `min(h, hₖ₋₁)` selects the larger |step| and |dt|
    # can grow on rejection instead of shrinking until the error test passes.
    prob = ODEProblem(fiip, [1.0], (1.0, 0.0), 1.0)
    for (alg, set_estm1) in (
            (FBDF(), (integ, v) -> (integ.cache.terkm1 = v)),
            (QNDF(), (integ, v) -> (integ.cache.EEst1 = v)),
        )
        integ = init(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        step!(integ)
        c = integ.cache
        c.order = 3
        c.consfailcnt = 5
        integ.dt = -1.0e-3
        OrdinaryDiffEqCore.set_EEst!(integ, 1.5)
        set_estm1(integ, 1.0e-3) # Fₖ₋₁ ≈ 7.7; signed min() would grow |dt| ~3.85×
        OrdinaryDiffEqCore.step_reject_controller!(integ, integ.alg)
        @test abs(integ.dt) <= 5.0e-4 # h was halved (cf > 1); |dt| must not exceed it
        @test c.order == 2
    end
end

# Regression test for the backward-in-time step rejection path: for tdir < 0
# (e.g. adjoint solves), `bdf_step_reject_controller!` must still shrink |dt|.
# Previously `min(h, hₖ₋₁)`/`hₖ₋₁ > hₖ` compared signed (negative) step sizes,
# picking the *larger* magnitude and letting dt hit a fixed point where
# EEst > 1 forever, hanging the solver at maxiters.
@testset "BDF step rejection shrinks |dt| for backward integration" begin
    for (t0, t1) in ((0.0, 1.0), (1.0, 0.0))
        prob = ODEProblem((u, p, t) -> -u, 1.0, (t0, t1))
        integ = init(prob, FBDF(); abstol = 1.0e-8, reltol = 1.0e-8)
        integ.cache.consfailcnt = 4
        integ.cache.order = 3
        OrdinaryDiffEqCore.set_EEst!(integ, 2.0)
        # small k-1 estimate makes the lower-order candidate much larger
        integ.cache.terkm1 = 0.01
        dt0 = integ.dt
        OrdinaryDiffEqCore.step_reject_controller!(integ, FBDF())
        @test signbit(integ.dt) == signbit(dt0)
        @test abs(integ.dt) <= abs(dt0) / 2
    end
end

@testset "QNDF2 step size control does not thrash (#4332)" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0))
    for alg in (QNDF2(), QBDF2())
        sol = solve(prob, alg, reltol = 1.0e-8, abstol = 1.0e-8)
        @test sol.retcode == ReturnCode.Success
        @test sol.stats.nreject < sol.stats.naccept / 10
        @test abs(sol.u[end] - exp(-10.0)) < 1.0e-6
    end
end

@testset "DFBDF first step: O(h²) BDF1 error estimate (#4791)" begin
    λ = 8.6e8
    u0 = [1.0e-3, 1.0e12]
    du0 = [-λ * u0[1], 1.0e-3 * u0[2]]
    fiip = (out, du, u, p, t) -> (out[1] = du[1] + λ * u[1]; out[2] = du[2] - 1.0e-3 * u[2]; nothing)
    foop = (du, u, p, t) -> [du[1] + λ * u[1], du[2] - 1.0e-3 * u[2]]
    tol = 1.0e-8
    for f in (fiip, foop)
        # Forced first step in the asymptotic regime: accepted with the analytic
        # error h²λ²u/2 (only u[1] contributes; RMS norm over two components).
        h = 1.8e-12
        prob = DAEProblem(f, du0, u0, (0.0, 1000.0); differential_vars = [true, true])
        integ = init(prob, DFBDF(); dt = h, abstol = tol, reltol = tol)
        step!(integ)
        @test integ.t - integ.tprev == h
        @test integ.stats.nreject == 0
        @test integ.u[1] ≈ u0[1] / (1 + h * λ) rtol = 1.0e-10
        EEst_exact = (h * λ)^2 * u0[1] / 2 / (tol + tol * u0[1]) / sqrt(2)
        @test OrdinaryDiffEqCore.get_EEst(integ) ≈ EEst_exact rtol = 0.01

        # At t0 = 14400 the first step must resolve the transient, and steps of a
        # few eps(t0) must not be rejected for using the uncommitted dt.
        t0 = 14400.0
        prob = DAEProblem(f, du0, u0, (t0, t0 + 1000.0); differential_vars = [true, true])
        for dt in (nothing, 1.0e-11)
            kw = dt === nothing ? (;) : (; dt)
            sol = solve(prob, DFBDF(); abstol = tol, reltol = tol, kw...)
            @test sol.retcode == ReturnCode.Success
            @test sol.u[end][2] ≈ u0[2] * exp(1.0) rtol = 1.0e-7
            # Dense output in the transient: an undetected O(h) first step is off by ~1e-4.
            for s in (1.0e-10, 1.0e-9, 1.0e-8)
                @test abs(sol(t0 + s)[1] - u0[1] * exp(-λ * s)) < 100 * tol
            end
        end
    end
end

@testset "DFBDF start without a consistent du0 (#4791)" begin
    function rober(out, du, u, p, t)
        out[1] = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
        out[2] = 0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
        out[3] = u[1] + u[2] + u[3] - 1.0
        return nothing
    end
    dv = [true, true, false]
    t0 = 1.0e4
    tspan = (t0, t0 + 100.0)
    ref = solve(
        DAEProblem(rober, [-0.04, 0.04, 0.0], [1.0, 0.0, 0.0], tspan; differential_vars = dv),
        DFBDF(); abstol = 1.0e-10, reltol = 1.0e-10
    )
    # The residual does not constrain the algebraic du0[3]; an inconsistent algebraic
    # u0[3] fails the residual check and keeps the old first step.
    for (du0, u0, init) in (
            ([-0.04, 0.04, -1.0e6], [1.0, 0.0, 0.0], CheckInit()),
            ([-0.04, 0.04, 0.0], [1.0, 0.0, 0.1], NoInit()),
        )
        prob = DAEProblem(rober, du0, u0, tspan; differential_vars = dv)
        sol = solve(prob, DFBDF(); abstol = 1.0e-8, reltol = 1.0e-8, initializealg = init)
        @test sol.retcode == ReturnCode.Success
        @test sol.u[end][1] ≈ ref.u[end][1] rtol = 1.0e-5
        @test sum(sol.u[end]) ≈ 1 atol = 1.0e-8
    end
end

@testset "DFBDF restart after a derivative discontinuity is a BDF1 cold start (#4791)" begin
    λ = 8.6e8
    u0 = [1.0e-3, 1.0e12]
    du0 = [-λ * u0[1], 1.0e-3 * u0[2]]
    fiip = (out, du, u, p, t) -> (out[1] = du[1] + λ * u[1]; out[2] = du[2] - 1.0e-3 * u[2]; nothing)
    foop = (du, u, p, t) -> [du[1] + λ * u[1], du[2] - 1.0e-3 * u[2]]
    tol = 1.0e-8
    prob(f, t0) = DAEProblem(f, du0, u0, (t0, t0 + 1000.0); differential_vars = [true, true])
    for f in (fiip, foop), t0 in (0.0, 14400.0), npre in (1, 4, 40), via_callback in (false, true)
        fire = Ref(false)
        cb = DiscreteCallback((u, t, integ) -> fire[], integ -> (fire[] = false; integ.u .= integ.u))
        integ = init(
            prob(f, t0), DFBDF(); abstol = tol, reltol = tol, callback = cb,
            initializealg = BrownFullBasicInit()
        )
        for _ in 1:npre
            step!(integ)
        end
        npre == 40 && @test integ.cache.order >= 2
        if via_callback
            fire[] = true
            step!(integ)
            @test !fire[]
        else
            derivative_discontinuity!(integ, true)
        end
        uprev = copy(integ.u)
        set_proposed_dt!(integ, max(1.0e-3 / λ, 2 * eps(integ.t)))
        nrej = integ.stats.nreject
        step!(integ)
        h = integ.t - integ.tprev
        @test integ.stats.nreject == nrej
        @test integ.u[1] ≈ uprev[1] / (1 + h * λ) rtol = 1.0e-8
        w = tol + tol * abs(uprev[1])
        @test OrdinaryDiffEqCore.get_EEst(integ) ≈ (h * λ)^2 * uprev[1] / 2 / w / sqrt(2) rtol = 0.01
    end
end
