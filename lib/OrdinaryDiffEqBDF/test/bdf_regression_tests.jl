using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, ForwardDiff, Test
using OrdinaryDiffEqCore: DEVerbosity
import OrdinaryDiffEqCore.SciMLLogging as SciMLLogging
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, NLNewton

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

@testset "MOOSE234 reinit" begin
    foop_m = (u, p, t) -> p * u
    proboop_m = ODEProblem(foop_m, ones(2), (0.0, 10.0), -1.0)

    fiip_m = (du, u, p, t) -> du .= p .* u
    probiip_m = ODEProblem(fiip_m, ones(2), (0.0, 10.0), -1.0)

    for prob in [proboop_m, probiip_m]
        integ = init(prob, MOOSE234())
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.sol.t[end] == 10.0

        reinit!(integ, prob.u0)
        solve!(integ)
        @test integ.sol.retcode == ReturnCode.Success
        @test integ.sol.t[end] == 10.0
    end
end

@testset "MOOSE234 startup order progression" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 5.0))
    integ = init(prob, MOOSE234(), dt = 0.01)

    @test integ.cache.order == 2
    @test integ.cache.iters_from_event == 0

    step!(integ)
    @test integ.cache.iters_from_event >= 1
    @test 2 <= integ.cache.order <= 4

    step!(integ)
    @test integ.cache.iters_from_event >= 2

    for _ in 1:10
        step!(integ)
    end
    @test integ.cache.iters_from_event >= 4
end

@testset "MOOSE234 adaptive order switching" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 10.0))
    integ = init(prob, MOOSE234(), dt = 0.01)

    orders_seen = Set{Int}()
    for _ in 1:200
        step!(integ)
        push!(orders_seen, integ.cache.order)
        integ.sol.retcode == ReturnCode.Failure && break
    end

    @test length(orders_seen) >= 2
    @test all(o -> 2 <= o <= 4, orders_seen)
end

@testset "MOOSE234 Van der Pol (stiff)" begin
    function vdp!(du, u, p, t)
        μ = p[1]
        du[1] = u[2]
        du[2] = μ * (1 - u[1]^2) * u[2] - u[1]
    end

    prob_mod = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3), [100.0])
    sol_mod = solve(prob_mod, MOOSE234(), abstol = 1e-6, reltol = 1e-6)
    @test sol_mod.retcode == ReturnCode.Success
    @test sol_mod.t[end] == 6.3

    prob_stiff = ODEProblem(vdp!, [2.0, 0.0], (0.0, 6.3), [1000.0])
    sol_stiff = solve(prob_stiff, MOOSE234(), abstol = 1e-4, reltol = 1e-4)
    @test sol_stiff.retcode == ReturnCode.Success
    @test sol_stiff.t[end] == 6.3
end

@testset "MOOSE234 ROBER stiff system" begin
    function rober!(du, u, p, t)
        y1, y2, y3 = u
        du[1] = -0.04 * y1 + 1e4 * y2 * y3
        du[2] = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2^2
        du[3] = 3e7 * y2^2
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1e5))
    sol = solve(prob, MOOSE234(), abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
    @test sol.t[end] == 1e5
    @test isapprox(sum(sol.u[end]), 1.0, atol = 1e-6)
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
