using OrdinaryDiffEqBDF, OrdinaryDiffEqRosenbrock, ODEProblemLibrary
using SciMLBase: DAEProblem, ODEProblem, ODEFunction, successful_retcode, remake,
    DiscreteCallback
using OrdinaryDiffEqNonlinearSolve: NLNewton, NonlinearSolveAlg
using LinearAlgebra, Test

@testset "NordsieckBDF: adaptive accuracy" begin
    for (nm, prob) in (
            ("out-of-place", ODEProblemLibrary.prob_ode_linear),
            ("in-place", ODEProblemLibrary.prob_ode_2Dlinear),
        )
        @testset "$nm" begin
            exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[2])
            prev = Inf
            for tol in (1.0e-6, 1.0e-8, 1.0e-10)
                sol = solve(prob, NordsieckBDF(), abstol = tol, reltol = tol)
                @test successful_retcode(sol)
                err = norm(sol.u[end] .- exact) / norm(exact)
                @test err < 100tol
                # tightening the tolerance must not make the answer worse
                @test err <= max(prev, 10eps())
                prev = err
            end
        end
    end
end

@testset "NordsieckBDF: max_order variants" begin
    prob = ODEProblemLibrary.prob_ode_2Dlinear
    exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[2])
    for mo in 1:5
        sol = solve(
            prob, NordsieckBDF(max_order = Val(mo)), abstol = 1.0e-8, reltol = 1.0e-8
        )
        @test successful_retcode(sol)
        # order 1 is legitimately much less accurate at a given tolerance
        @test norm(sol.u[end] .- exact) / norm(exact) < (mo == 1 ? 1.0e-3 : 1.0e-5)
        @test sol.stats.naccept > 0
    end
end

@testset "NordsieckBDF: stiff problems" begin
    probs = (
        ("ROBER", remake(ODEProblemLibrary.prob_ode_rober, tspan = (0.0, 1.0e5))),
        ("HIRES", ODEProblemLibrary.prob_ode_hires),
        ("OREGO", ODEProblemLibrary.prob_ode_orego),
        ("POLLU", ODEProblemLibrary.prob_ode_pollution),
        ("VDPOL", ODEProblemLibrary.prob_ode_vanderpol_stiff),
    )
    for (nm, prob) in probs
        @testset "$nm" begin
            ref = solve(prob, Rodas5P(), abstol = 1.0e-14, reltol = 1.0e-14)
            rv = ref(prob.tspan[2])
            sol = solve(
                prob, NordsieckBDF(), abstol = 1.0e-10, reltol = 1.0e-8,
                save_everystep = false
            )
            @test successful_retcode(sol)
            @test norm(sol.u[end] .- rv) / norm(rv) < 1.0e-4
        end
    end
end

@testset "NordsieckBDF: dense output is the Nordsieck polynomial" begin
    prob = ODEProblemLibrary.prob_ode_hires
    ref = solve(prob, Rodas5P(), abstol = 1.0e-14, reltol = 1.0e-14)
    sol = solve(prob, NordsieckBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
    worst = 0.0
    for t in range(prob.tspan[1], prob.tspan[2], length = 53)
        worst = max(worst, norm(sol(t) .- ref(t)) / max(norm(ref(t)), 1.0e-10))
    end
    @test worst < 1.0e-4
    # interpolating at the step points must reproduce the stored solution
    @test sol(sol.t[end]) ≈ sol.u[end]
    @test sol(sol.t[1]) ≈ sol.u[1]
end

@testset "NordsieckBDF: mass matrix" begin
    # M u' = A u with M singular-free (index-0), analytic solution via M\A
    M = [2.0 0.0; 0.0 4.0]
    A = [-1.0 0.5; 0.25 -2.0]
    f = ODEFunction((du, u, p, t) -> mul!(du, A, u), mass_matrix = M)
    u0 = [1.0, 2.0]
    prob = ODEProblem(f, u0, (0.0, 1.0))
    sol = solve(prob, NordsieckBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test successful_retcode(sol)
    ref = solve(prob, Rodas5P(), abstol = 1.0e-13, reltol = 1.0e-13)
    @test norm(sol.u[end] .- ref.u[end]) / norm(ref.u[end]) < 1.0e-6
end

@testset "NordsieckBDF: callback restarts the history" begin
    # a discontinuity forces the Nordsieck array back to order 1
    f = (du, u, p, t) -> (du[1] = -u[1]; nothing)
    prob = ODEProblem(f, [1.0], (0.0, 2.0))
    cb = DiscreteCallback((u, t, integ) -> t == 1.0, integ -> (integ.u[1] += 1.0))
    sol = solve(
        prob, NordsieckBDF(), abstol = 1.0e-10, reltol = 1.0e-10,
        callback = cb, tstops = [1.0]
    )
    @test successful_retcode(sol)
    # analytic: exp(-t) up to t=1, then (exp(-1)+1)exp(-(t-1))
    want = (exp(-1.0) + 1.0) * exp(-1.0)
    @test isapprox(sol.u[end][1], want, rtol = 1.0e-5)
end

@testset "DNordsieckBDF: index-1 DAE" begin
    # u1' = -u1,  0 = u2 - u1^2  =>  u1 = e^{-t}, u2 = e^{-2t}
    function dae_res!(res, du, u, p, t)
        res[1] = du[1] + u[1]
        res[2] = u[2] - u[1]^2
        return nothing
    end
    dae_res_oop(du, u, p, t) = [du[1] + u[1], u[2] - u[1]^2]
    u0 = [1.0, 1.0]
    du0 = [-1.0, -2.0]
    dvars = [true, false]
    exact = [exp(-2.0), exp(-4.0)]

    prob_iip = DAEProblem(dae_res!, du0, u0, (0.0, 2.0), differential_vars = dvars)
    prob_oop = DAEProblem{false}(
        dae_res_oop, du0, u0, (0.0, 2.0), differential_vars = dvars
    )
    for (nm, prob) in (("in-place", prob_iip), ("out-of-place", prob_oop))
        @testset "$nm" begin
            prev = Inf
            for tol in (1.0e-6, 1.0e-8, 1.0e-10)
                sol = solve(prob, DNordsieckBDF(), abstol = tol, reltol = tol)
                @test successful_retcode(sol)
                err = norm(sol.u[end] .- exact) / norm(exact)
                @test err < 100tol
                @test err <= max(prev, 10eps())
                prev = err
                # the algebraic constraint is enforced to the solver tolerance
                @test isapprox(
                    sol.u[end][2], sol.u[end][1]^2, rtol = max(1.0e-7, 1000tol)
                )
            end
        end
    end
end

@testset "DNordsieckBDF: Robertson DAE" begin
    function rober_dae!(out, du, u, p, t)
        out[1] = -0.04u[1] + 1.0e4 * u[2] * u[3] - du[1]
        out[2] = +0.04u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3] - du[2]
        out[3] = u[1] + u[2] + u[3] - 1.0
        return nothing
    end
    prob = DAEProblem(
        rober_dae!, [-0.04, 0.04, 0.0], [1.0, 0.0, 0.0], (0.0, 1.0e5),
        differential_vars = [true, true, false]
    )
    sol = solve(prob, DNordsieckBDF(), abstol = 1.0e-8, reltol = 1.0e-8)
    @test successful_retcode(sol)
    # mass conservation is the algebraic constraint
    @test isapprox(sum(sol.u[end]), 1.0, atol = 1.0e-10)
    ref = solve(prob, DFBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test norm(sol.u[end] .- ref.u[end]) / norm(ref.u[end]) < 1.0e-4
end

@testset "NordsieckBDF: loose corrector keeps the step size" begin
    # The point of the Nordsieck representation: loosening the nonlinear solve to a
    # fraction of the local error budget must not blow up the step count, which is
    # what happens when the predictor is rebuilt from stored history.
    prob = ODEProblemLibrary.prob_ode_hires
    tight = solve(
        prob, NordsieckBDF(nlsolve = NLNewton(κ = 1 // 100)),
        abstol = 1.0e-8, reltol = 1.0e-6, save_everystep = false
    )
    loose = solve(
        prob, NordsieckBDF(nlsolve = NLNewton(κ = 1 // 10)),
        abstol = 1.0e-8, reltol = 1.0e-6, save_everystep = false
    )
    @test successful_retcode(tight) && successful_retcode(loose)
    ntight = tight.stats.naccept + tight.stats.nreject
    nloose = loose.stats.naccept + loose.stats.nreject
    @test nloose < 3 * ntight
end

@testset "NordsieckBDF: NonlinearSolveAlg backend" begin
    # The pluggable NonlinearSolve.jl corrector must give the same answer as the
    # built-in Newton, and `has_special_newton_error` is honoured on both paths so
    # `κ` keeps its NLSCOEF meaning either way.
    for (nm, prob) in (
            ("out-of-place", ODEProblemLibrary.prob_ode_linear),
            ("in-place", ODEProblemLibrary.prob_ode_2Dlinear),
        )
        @testset "$nm" begin
            exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[2])
            ref = solve(prob, NordsieckBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
            sol = solve(
                prob, NordsieckBDF(nlsolve = NonlinearSolveAlg()),
                abstol = 1.0e-10, reltol = 1.0e-10
            )
            @test successful_retcode(sol)
            @test norm(sol.u[end] .- exact) / norm(exact) < 1.0e-6
            @test isapprox(sol.u[end], ref.u[end], rtol = 1.0e-6)
        end
    end

    # stiff problem, to exercise the W-reuse path
    prob = ODEProblemLibrary.prob_ode_hires
    sol = solve(
        prob, NordsieckBDF(nlsolve = NonlinearSolveAlg()),
        abstol = 1.0e-10, reltol = 1.0e-8, save_everystep = false
    )
    @test successful_retcode(sol)
    ref = solve(prob, Rodas5P(), abstol = 1.0e-14, reltol = 1.0e-14)
    @test norm(sol.u[end] .- ref(prob.tspan[2])) / norm(ref(prob.tspan[2])) < 1.0e-4
end

@testset "NordsieckBDF: stald and repeated-failure order reduction" begin
    # STALD is off by default (as in CVODE). Enabling it must not change results on
    # problems that are not stability-limited, and must not break anything.
    for prob in (
            ODEProblemLibrary.prob_ode_hires,
            remake(ODEProblemLibrary.prob_ode_rober, tspan = (0.0, 1.0e5)),
        )
        off = solve(
            prob, NordsieckBDF(), abstol = 1.0e-8, reltol = 1.0e-6,
            save_everystep = false
        )
        on = solve(
            prob, NordsieckBDF(stald = true), abstol = 1.0e-8, reltol = 1.0e-6,
            save_everystep = false
        )
        @test successful_retcode(off) && successful_retcode(on)
        @test isapprox(off.u[end], on.u[end], rtol = 1.0e-8)
    end
    @test successful_retcode(
        solve(
            DAEProblem(
                (res, du, u, p, t) -> (res[1] = du[1] + u[1]; nothing),
                [-1.0], [1.0], (0.0, 1.0), differential_vars = [true]
            ),
            DNordsieckBDF(stald = true), abstol = 1.0e-8, reltol = 1.0e-8
        )
    )

    # A jump in the RHS forces corrector failures; the solver must recover rather
    # than only shrinking dt forever.
    f = (du, u, p, t) -> (du[1] = t < 1.0 ? -u[1] : -100 * u[1] + 99; nothing)
    prob = ODEProblem(f, [1.0], (0.0, 2.0))
    sol = solve(prob, NordsieckBDF(), abstol = 1.0e-10, reltol = 1.0e-10)
    @test successful_retcode(sol)
    @test isapprox(sol.u[end][1], 0.99, atol = 1.0e-3)
end
