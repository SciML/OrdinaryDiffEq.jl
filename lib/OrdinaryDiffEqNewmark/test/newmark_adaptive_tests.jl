using OrdinaryDiffEqNewmark, Test, RecursiveArrayTools
using SciMLBase: ReturnCode

# Harmonic oscillator with analytic solution
f1_harmonic_iip!(dv, v, u, p, t) = (dv .= -u)
f2_harmonic_iip!(du, v, u, p, t) = (du .= v)
f1_harmonic_oop(v, u, p, t) = -u
f2_harmonic_oop(v, u, p, t) = v
function harmonic_analytic(y0, p, t)
    v0, u0 = y0.x
    return ArrayPartition(-u0 * sin(t) + v0 * cos(t), u0 * cos(t) + v0 * sin(t))
end

function harmonic_probs(N, T = 5.0)
    v0 = ones(N)
    u0 = zeros(N)
    iip = DynamicalODEProblem(
        DynamicalODEFunction(
            f1_harmonic_iip!, f2_harmonic_iip!; analytic = harmonic_analytic
        ),
        v0, u0, (0.0, T)
    )
    oop = DynamicalODEProblem(
        DynamicalODEFunction(
            f1_harmonic_oop, f2_harmonic_oop; analytic = harmonic_analytic
        ),
        v0, u0, (0.0, T)
    )
    return iip, oop
end

function harmonic_maxerr(sol)
    m = 0.0
    for (i, t) in enumerate(sol.t)
        ex = harmonic_analytic(sol.prob.u0, nothing, t)
        m = max(m, maximum(abs.(sol.u[i] .- ex)))
    end
    return m
end

# Nonlinear pendulum
f1_pendulum_iip!(dv, v, u, p, t) = (dv .= -9.81 .* sin.(u))
f2_pendulum_iip!(du, v, u, p, t) = (du .= v)
pendulum_prob() = DynamicalODEProblem(
    DynamicalODEFunction(f1_pendulum_iip!, f2_pendulum_iip!),
    [0.0], [0.5], (0.0, 5.0)
)

@testset "Adaptive solves with default options" begin
    for N in (1, 2, 3)
        prob_iip, prob_oop = harmonic_probs(N)
        for prob in (prob_iip, prob_oop)
            for alg in (NewmarkBeta(), GeneralizedAlpha(rho_inf = 0.9))
                sol = solve(prob, alg, dt = 0.1)
                @test sol.retcode == ReturnCode.Success
                # default tolerances are reltol = 1e-3, so a few 1e-2 of
                # accumulated error over five periods is in spec
                @test harmonic_maxerr(sol) < 1.0e-1
            end
        end
    end
end

@testset "Error tracks tolerance monotonically" begin
    prob_iip, prob_oop = harmonic_probs(2)
    for prob in (prob_iip, prob_oop)
        errs = Float64[]
        nsteps = Int[]
        for tol in (1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12)
            sol = solve(prob, NewmarkBeta(), dt = 0.1, abstol = tol, reltol = tol)
            @test sol.retcode == ReturnCode.Success
            push!(errs, harmonic_maxerr(sol))
            push!(nsteps, length(sol.t) - 1)
            # The per-step error target is tol, so the global error on this
            # non-dissipative problem accumulates to O(nsteps * tol); healthy
            # baseline is Q ≈ 1-3.
            @test errs[end] <= 10 * nsteps[end] * tol
        end
        @test issorted(errs, rev = true)
        # While the position channel binds, dt ~ tol^(1/3) and 100x tighter
        # tolerance costs 100^(1/3) ≈ 4.6x the steps; once the clamped
        # velocity surrogate binds (γ = 1/2 is inside the velocity band),
        # dt ~ tol^(1/2) and the factor approaches 100^(1/2) = 10.
        for i in 1:(length(nsteps) - 1)
            @test 3 < nsteps[i + 1] / nsteps[i] < 12
        end
        # global error ~ (mean dt)^2: adaptive attains the method order
        slope = (log(errs[1]) - log(errs[end])) / (log(nsteps[end]) - log(nsteps[1]))
        @test 1.7 < slope < 2.3
    end
end

@testset "Uninformative estimator band is rejected" begin
    prob_iip, prob_oop = harmonic_probs(2)
    for β in (
            1 / 6, nextfloat(1 / 6), prevfloat(1 / 6),
            1 / 6 + 1.0e-6, 1 / 6 - 1.0e-6, 1 / 6 + 1.0e-8, 1 / 6 - 1.0e-8,
            1 / 6 + 1.0e-10, 1 / 6 - 1.0e-10, 0.1667,
        )
        for prob in (prob_iip, prob_oop)
            @test_throws ArgumentError solve(prob, NewmarkBeta(β, 0.5), dt = 0.1)
            # the parameters stay usable with fixed steps
            sol = solve(prob, NewmarkBeta(β, 0.5), dt = 0.1, adaptive = false)
            @test sol.retcode == ReturnCode.Success
        end
    end
    # γ within its band around 1/2 while β is inside its own band: rejected
    @test_throws ArgumentError solve(
        prob_iip, NewmarkBeta(1 / 6, 0.5 + 1.0e-3), dt = 0.1
    )
    # a healthy velocity channel rescues a vanishing position channel
    for prob in (prob_iip, prob_oop)
        sol = solve(prob, NewmarkBeta(1 / 6, 0.7), dt = 0.1)
        @test sol.retcode == ReturnCode.Success
        @test harmonic_maxerr(sol) < 1.0e-1
    end
    # generalized-α: the loci move to β = 1/(6κ), γ = 1/(2κ), κ = (1-αf)/(1-αm)
    κ = (1 - 0.2) / (1 - 0.1)
    @test_throws ArgumentError solve(
        prob_iip, GeneralizedAlpha(0.1, 0.2, 1 / (6κ), 1 / (2κ)), dt = 0.1
    )
    @test_throws ArgumentError solve(
        prob_oop, GeneralizedAlpha(0.1, 0.2, 1 / (6κ), 1 / (2κ)), dt = 0.1
    )
end

@testset "Degenerate estimate cannot grow dt" begin
    # Constant acceleration: aₙ₊₁ - aₙ = 0 on every step, so the estimate is
    # exactly zero; the controller must hold the step size instead of walking
    # it out to dtmax.
    g1!(dv, v, u, p, t) = (dv .= -9.81)
    g2!(du, v, u, p, t) = (du .= v)
    gprob = DynamicalODEProblem(
        DynamicalODEFunction(g1!, g2!), [1.0], [0.0], (0.0, 10.0)
    )
    sol = solve(gprob, NewmarkBeta(), dt = 0.5, dtmax = 1.0e6)
    @test sol.retcode == ReturnCode.Success
    @test maximum(diff(sol.t)) <= 0.5 + 1.0e-12
    @test abs(sol.u[end].x[2][1] - (10.0 - 9.81 * 50.0)) < 1.0e-8
    # same defence on the out-of-place estimate path
    g1(v, u, p, t) = [-9.81]
    g2(v, u, p, t) = v
    gprob_oop = DynamicalODEProblem(
        DynamicalODEFunction(g1, g2), [1.0], [0.0], (0.0, 10.0)
    )
    sol = solve(gprob_oop, NewmarkBeta(), dt = 0.5, dtmax = 1.0e6)
    @test sol.retcode == ReturnCode.Success
    @test maximum(diff(sol.t)) <= 0.5 + 1.0e-12
    @test abs(sol.u[end].x[2][1] - (10.0 - 9.81 * 50.0)) < 1.0e-8
end

@testset "Velocity surrogate controls high-frequency error" begin
    # u'' = -ω²u with unit velocity amplitude: the position tolerance scale is
    # 1/ω smaller than the velocity scale, so a position-only estimate lets
    # the velocity error through by a factor ~ω. Before the clamped velocity
    # channel this case measured Q = verr/(nsteps*tol) ≈ 119 with a Success
    # return; the surrogate holds it at ≈ 0.4.
    ω = 100.0
    s1!(dv, v, u, p, t) = (dv .= -ω^2 .* u)
    s2!(du, v, u, p, t) = (du .= v)
    s1(v, u, p, t) = -ω^2 .* u
    s2(v, u, p, t) = v
    tol = 1.0e-6
    for (fns, alg) in (
            ((s1!, s2!), NewmarkBeta()),
            ((s1, s2), NewmarkBeta()),
            ((s1!, s2!), GeneralizedAlpha(rho_inf = 0.9)),
        )
        prob = DynamicalODEProblem(
            DynamicalODEFunction(fns...), [1.0], [0.0], (0.0, 5.0)
        )
        sol = solve(
            prob, alg, dt = 0.01, abstol = tol, reltol = tol, maxiters = 10^7
        )
        @test sol.retcode == ReturnCode.Success
        verr = maximum(
            abs(sol.u[i].x[1][1] - cos(ω * t)) for (i, t) in enumerate(sol.t)
        )
        nsteps = length(sol.t) - 1
        @test verr <= 20 * nsteps * tol
    end
end

# Pins the position channel itself: with it zeroed, this configuration takes
# 19 steps at Q ≈ 5.6 instead of 34 at Q ≈ 0.7, so either assertion fails.
@testset "Position channel binds the step size" begin
    prob = DynamicalODEProblem(
        DynamicalODEFunction(
            (dv, v, u, p, t) -> (dv .= -u), (du, v, u, p, t) -> (du .= v)
        ), [0.0], [1.0], (0.0, 5.0)
    )
    tol = 1.0e-3
    sol = solve(prob, NewmarkBeta(0.5, 0.5), dt = 0.1, abstol = tol, reltol = tol)
    @test sol.retcode == ReturnCode.Success
    @test sol.stats.naccept >= 30
    nsteps = length(sol.t) - 1
    maxerr = maximum(
        abs(sol.u[i].x[2][1] - cos(t)) for (i, t) in enumerate(sol.t)
    )
    @test maxerr <= 2 * nsteps * tol
end

@testset "Out-of-place saveat interpolation" begin
    prob = DynamicalODEProblem(
        DynamicalODEFunction(f1_harmonic_oop, f2_harmonic_oop),
        zeros(2), ones(2), (0.0, 1.0)
    )
    # GeneralizedAlpha(rho_inf = 1) runs the same update arithmetic as
    # NewmarkBeta() but through the GeneralizedAlphaConstantCache step, whose
    # uprev-aliasing behaviour is otherwise untested.
    for alg in (NewmarkBeta(), GeneralizedAlpha(rho_inf = 1.0))
        sol = solve(prob, alg, dt = 0.1, saveat = 0.25, adaptive = false)
        @test sol.retcode == ReturnCode.Success
        for (i, t) in enumerate(sol.t)
            @test maximum(abs.(sol.u[i].x[1] .- (-sin(t)))) < 1.0e-3
            @test maximum(abs.(sol.u[i].x[2] .- cos(t))) < 1.0e-3
        end
    end
end

@testset "Fixed-step results are unchanged" begin
    # Damped oscillator, dt = 0.1, NewmarkBeta(): the final state is pinned
    # bit-for-bit to the values the solver produced before adaptive support.
    d1!(dv, v, u, p, t) = (dv .= -u .- 0.3 .* v)
    d2!(du, v, u, p, t) = (du .= v)
    d1(v, u, p, t) = -u .- 0.3 .* v
    d2(v, u, p, t) = v
    v0 = [1.0, 0.3]
    u0 = [0.0, 0.7]
    prob_iip = DynamicalODEProblem(DynamicalODEFunction(d1!, d2!), v0, u0, (0.0, 3.0))
    prob_oop = DynamicalODEProblem(DynamicalODEFunction(d1, d2), v0, u0, (0.0, 3.0))
    anchor = reinterpret.(
        Float64,
        [0xbfe4a7055bdc734a, 0xbfd181aebd616175, 0x3fbd3aaf9181c633, 0xbfd92fc96ff9344d]
    )
    sol_iip = solve(prob_iip, NewmarkBeta(), dt = 0.1, adaptive = false)
    sol_oop = solve(prob_oop, NewmarkBeta(), dt = 0.1, adaptive = false)
    @test collect(sol_iip.u[end]) == anchor
    @test collect(sol_oop.u[end]) == anchor
    # in-place and out-of-place steps run the same arithmetic
    @test all(sol_iip.u[i] == sol_oop.u[i] for i in eachindex(sol_iip.u))
end

@testset "Parameter sweep stays within tolerance" begin
    prob_iip, prob_oop = harmonic_probs(2)
    in_band(β, γ) = abs(β - 1 / 6) < 1 / 24 && abs(γ - 1 / 2) < 1 / 100
    n_success = 0
    for β in (0.05, 1 / 6, 0.25, 1 / 3, 0.5)
        for γ in (0.3, 0.5, 0.7, 1.0)
            for tol in (1.0e-6, 1.0e-10)
                for prob in (prob_iip, prob_oop)
                    alg = NewmarkBeta(β, γ)
                    if in_band(β, γ)
                        @test_throws ArgumentError solve(
                            prob, alg, dt = 0.1, abstol = tol, reltol = tol
                        )
                        continue
                    end
                    sol = solve(
                        prob, alg, dt = 0.1, abstol = tol, reltol = tol,
                        maxiters = 10^7
                    )
                    if sol.retcode == ReturnCode.Success
                        n_success += 1
                        nsteps = length(sol.t) - 1
                        @test harmonic_maxerr(sol) <= 20 * nsteps * tol
                    end
                end
            end
        end
    end
    # The accuracy assertion above only runs on Success, so a regression that
    # turns solves into failures would silently hollow the sweep out; pin the
    # Success count (all 76 non-rejected combinations succeed today).
    @test n_success == 76
    for ρ in (0.4, 0.7, 1.0)
        for tol in (1.0e-6, 1.0e-8)
            sol = solve(
                prob_iip, GeneralizedAlpha(rho_inf = ρ), dt = 0.1,
                abstol = tol, reltol = tol, maxiters = 10^7
            )
            @test sol.retcode == ReturnCode.Success
            @test harmonic_maxerr(sol) <= 20 * (length(sol.t) - 1) * tol
        end
    end
    # ρ∞ < 1/3 implies γ > 1 and has never been constructible
    @test_throws AssertionError GeneralizedAlpha(rho_inf = 0.2)

    # Nonlinear pendulum against a tight-tolerance reference. Q = err/(n*tol)
    # is problem-dependent: the harmonic-oscillator baseline is Q ≈ 0.2-3,
    # while this near-separatrix pendulum (u0 = 2.5 rad) measures Q ≈ 34-39 for the velocity-channel
    # configuration NewmarkBeta(1/6, 0.7), which is why only the default and
    # ρ∞ = 0.9 configurations are asserted against the 20x bound here.
    pprob = pendulum_prob()
    ref = solve(
        pprob, NewmarkBeta(), dt = 0.1, abstol = 1.0e-13, reltol = 1.0e-13,
        save_everystep = false, maxiters = 10^7
    )
    for tol in (1.0e-6, 1.0e-10)
        for alg in (NewmarkBeta(), GeneralizedAlpha(rho_inf = 0.9))
            sol = solve(
                pprob, alg, dt = 0.1, abstol = tol, reltol = tol,
                save_everystep = false, maxiters = 10^7
            )
            @test sol.retcode == ReturnCode.Success
            @test maximum(abs.(sol.u[end] .- ref.u[end])) <= 20 * sol.stats.naccept * tol
        end
    end
end
