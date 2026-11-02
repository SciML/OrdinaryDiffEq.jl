using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqRosenbrock, SciMLBase, Test

# u' = -p u has the exact solution u0 exp(-p (t - t0)). The decay is slow, so a run that
# continues from stale history stays visibly wrong instead of relaxing back onto the solution.
decay!(du, u, p, t) = (du[1] = -p[1] * u[1]; nothing)
decay(u, p, t) = [-p[1] * u[1]]
decay_dae!(r, du, u, p, t) = (r[1] = du[1] + p[1] * u[1]; nothing)
decay_dae(du, u, p, t) = [du[1] + p[1] * u[1]]

const TOLS = (; abstol = 1.0e-8, reltol = 1.0e-8)
exact_error(sol, u0, t0, tf, p) = maximum(
    t -> abs(sol(t)[1] - u0 * exp(-p * (t - t0))), range(t0, tf; length = 41)
)

function run_to_end_or_midway!(integ, stop)
    if stop == :end
        solve!(integ)
    else
        while integ.t < 5
            step!(integ)
        end
    end
    return integ
end

function test_matches_fresh(integ, fresh, u0, t0, tf, p, tol)
    @test integ.sol.retcode == ReturnCode.Success
    @test integ.sol.t == fresh.t
    @test integ.sol.u == fresh.u
    @test exact_error(integ.sol, u0, t0, tf, p) < tol
    return nothing
end

const ODE_ALGS = (
    FBDF => 1.0e-6, QNDF => 1.0e-6, QBDF => 1.0e-6, QNDF1 => 1.0e-3, QBDF1 => 1.0e-3,
    QNDF2 => 1.0e-4, QBDF2 => 1.0e-4, ABDF2 => 1.0e-3, NordsieckBDF => 1.0e-6,
)

@testset "reinit! restarts the multistep history of $alg" for (alg, tol) in ODE_ALGS
    @testset "inplace=$iip stop=$stop" for iip in (true, false), stop in (:end, :midway)
        f = iip ? decay! : decay
        prob = ODEProblem(f, [1.0], (0.0, 10.0), [1.0])
        fresh = solve(prob, alg(); TOLS...)
        integ = init(prob, alg(); TOLS...)

        run_to_end_or_midway!(integ, stop)
        reinit!(integ)
        solve!(integ)
        test_matches_fresh(integ, fresh, 1.0, 0.0, 10.0, 1.0, tol)

        reinit!(integ)
        solve!(integ)
        test_matches_fresh(integ, fresh, 1.0, 0.0, 10.0, 1.0, tol)

        new_state = solve(remake(prob; u0 = [3.0]), alg(); TOLS...)
        reinit!(integ, [3.0])
        solve!(integ)
        test_matches_fresh(integ, new_state, 3.0, 0.0, 10.0, 1.0, 3 * tol)

        new_span = solve(remake(prob; u0 = [2.0], tspan = (1.0, 7.0)), alg(); TOLS...)
        reinit!(integ, [2.0]; t0 = 1.0, tf = 7.0)
        solve!(integ)
        test_matches_fresh(integ, new_span, 2.0, 1.0, 7.0, 1.0, 2 * tol)

        new_p = solve(remake(prob; p = [2.0]), alg(); TOLS...)
        integ.p[1] = 2.0
        reinit!(integ; tf = 10.0)
        solve!(integ)
        test_matches_fresh(integ, new_p, 1.0, 0.0, 10.0, 2.0, tol)
    end
end

@testset "reinit! restarts the multistep history of $alg" for alg in (DFBDF, DNordsieckBDF)
    @testset "inplace=$iip" for iip in (true, false)
        F = iip ? decay_dae! : decay_dae
        prob = DAEProblem(F, [-1.0], [1.0], (0.0, 10.0), [1.0]; differential_vars = [true])
        fresh = solve(prob, alg(); TOLS...)
        integ = init(prob, alg(); TOLS...)
        solve!(integ)
        reinit!(integ)
        solve!(integ)
        test_matches_fresh(integ, fresh, 1.0, 0.0, 10.0, 1.0, 1.0e-6)
    end
end

# The switch into a composite member goes through `initialize!`, so the member that
# `reinit!` does not start on must not read the history of the previous solve either.
@testset "reinit! with a composite member entered at a switch: $alg" for (alg, tol) in (
        QNDF2 => 1.0e-4, QBDF2 => 1.0e-4, QNDF1 => 1.0e-3,
    )
    choice(integrator) = 2 <= integrator.t < 4 || 6 <= integrator.t < 8 ? 2 : 1
    composite() = OrdinaryDiffEqCore.CompositeAlgorithm((Rodas5P(), alg()), choice)
    prob = ODEProblem(decay!, [1.0], (0.0, 10.0), [1.0])
    fresh = solve(prob, composite(); TOLS...)
    @test length(unique(fresh.alg_choice)) == 2
    @test exact_error(fresh, 1.0, 0.0, 10.0, 1.0) < tol
    integ = init(prob, composite(); TOLS...)
    solve!(integ)
    reinit!(integ)
    solve!(integ)
    test_matches_fresh(integ, fresh, 1.0, 0.0, 10.0, 1.0, tol)
end
