using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, DiffEqBase, ForwardDiff, Test, LinearAlgebra
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

# choose_order! recomputes terkm2 while lowering; that value enters the loop
# condition on every drop from order ≥ 4. On u = exp(-t), terkm2 at order k
# approximates ‖h^(k-2) u^(k-2)‖ when t and dt match the FD stencil.
@testset "choose_order! terkm2 matches analytic h^(k-2) u^(k-2) on exp(-t)" begin
    choose_order! = OrdinaryDiffEqBDF.choose_order!
    calculate_residuals = DiffEqBase.calculate_residuals

    function residual_norm(integ, terk_tmp)
        atmp = calculate_residuals(
            terk_tmp, integ.uprev, integ.u,
            integ.opts.abstol, integ.opts.reltol,
            integ.opts.internalnorm, integ.t
        )
        return integ.opts.internalnorm(atmp, integ.t)
    end

    # After step!, restore the pre-advance call-site state used by choose_order!:
    # t is the step start and t + dt == ts_tmp[1].
    function restore_call_site!(integ)
        cache = integ.cache
        t_end = cache.ts_tmp[1]
        t_start = cache.ts_tmp[2]
        integ.dt = t_end - t_start
        integ.t = t_start
        return nothing
    end

    dae_exp = DAEProblem(
        (res, du, u, p, t) -> (res[1] = du[1] + u[1]),
        [-1.0], [1.0], (0.0, 2.0);
        differential_vars = [true]
    )

    for (prob, alg) in (
            (ODEProblem((u, p, t) -> -u, 1.0, (0.0, 2.0)), FBDF()),
            (ODEProblem((du, u, p, t) -> (du[1] = -u[1]), [1.0], (0.0, 2.0)), FBDF()),
            (dae_exp, DFBDF()),
        )
        integ = init(prob, alg; abstol = 1.0e-10, reltol = 0.0, dt = 1.0e-2)
        while integ.cache.order < 5 && integ.t < 1.0
            step!(integ)
        end
        @test integ.cache.order >= 5
        restore_call_site!(integ)
        @test integ.t + integ.dt ≈ integ.cache.ts_tmp[1] atol = 1.0e-14

        cache = integ.cache
        cache.order = 5
        cache.qwait = 1 # block order raise
        # Non-monotone terk chain forces lowering; from order 5 the returned
        # terk is the first recomputed terkm2, i.e. terkm2 at order 4.
        cache.terkm2 = 1.0
        cache.terkm1 = 2.0
        cache.terk = 3.0
        cache.terkp1 = 4.0

        # At order k = 4, terkm2 ≈ h^(k-2) u^(k-2) = h^2 u''. For u = exp(-t),
        # u'' = exp(-t), evaluated at the stencil point t + dt.
        k = 4
        m = k - 2
        tdt = integ.t + integ.dt
        analytic = (
            integ.u isa Number ? one(integ.u) :
                ones(eltype(integ.u), size(integ.u))
        ) *
            (abs(integ.dt)^m * exp(-tdt))
        analytic_res = residual_norm(integ, analytic)

        knew, terk = choose_order!(alg, integ, cache, Val(5))
        @test knew == 2
        @test terk ≈ analytic_res rtol = 1.0e-2
    end
end
