using OrdinaryDiffEqBDF, ForwardDiff, Test

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
    for MO in 1:5
        for prob in [proboop, probiip]
            sol = solve(prob, FBDF(max_order = Val{MO}()), verbose = false)
            @test sol.t[end] >= 700
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
