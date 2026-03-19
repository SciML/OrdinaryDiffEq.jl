using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, ForwardDiff, Test

foop = (u, p, t) -> u * p
proboop = ODEProblem(foop, ones(2), (0.0, 1000.0), 1.0)

fiip = (du, u, p, t) -> du .= u .* p
probiip = ODEProblem(fiip, ones(2), (0.0, 1000.0), 1.0)

@testset "FBDF reinit" begin
    for prob in [proboop, probiip]
        integ = init(prob, FBDF(), verbose = false) #suppress warning to clean up CI
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
        return sum(sol)
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
            sol = solve(prob, FBDF(max_order = Val{MO}()), verbose = false)
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
            dae_prob_mo, DFBDF(max_order = Val{MO}()),
            abstol = 1.0e-8, reltol = 1.0e-8, verbose = false
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
end # VERSION >= v"1.11"

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
            save_everystep = false
        )
        for _ in 1:10
            step!(integrator)
        end
        allocs = @allocated step!(integrator)
        @test allocs == 0
    end
end
