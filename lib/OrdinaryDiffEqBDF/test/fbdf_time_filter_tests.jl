using OrdinaryDiffEqBDF
using SciMLBase
using LinearAlgebra
using Test

exp_decay!(du, u, p, t) = (du[1] = -u[1]; nothing)
exp_decay(u, p, t) = -u

const IIP_DECAY = ODEProblem(exp_decay!, [1.0], (0.0, 10.0))
const OOP_DECAY = ODEProblem(exp_decay, 1.0, (0.0, 10.0))

# Largest deviation from exp(-t) over the saved solution, in units of the mixed
# tolerance the solve was run at, so it stays meaningful once exp(-t) drops below
# `abstol`.
function decay_error(sol, tol)
    return maximum(
        abs(only(u) - exp(-t)) / (tol + tol * exp(-t)) for (u, t) in zip(sol.u, sol.t)
    )
end

@testset "time_filter solves accurately" begin
    for (name, prob) in (("in-place", IIP_DECAY), ("out-of-place", OOP_DECAY))
        @testset "$name" begin
            for tol in (1.0e-6, 1.0e-8, 1.0e-10)
                sol = solve(
                    prob, FBDF(time_filter = true), abstol = tol, reltol = tol
                )
                @test SciMLBase.successful_retcode(sol)
                @test decay_error(sol, tol) < 100
            end
        end
    end
end

@testset "time_filter does not need more steps on a smooth problem" begin
    for prob in (IIP_DECAY, OOP_DECAY)
        for tol in (1.0e-6, 1.0e-8)
            filtered = solve(
                prob, FBDF(time_filter = true), abstol = tol, reltol = tol
            )
            plain = solve(prob, FBDF(), abstol = tol, reltol = tol)
            @test SciMLBase.successful_retcode(filtered)
            # The filter is only worth its extra rhs evaluations if it enlarges the
            # steps: it needs about 0.85x the steps here. Slack so that this is a
            # sanity check rather than a step-count regression test.
            @test length(filtered.t) <= 1.1 * length(plain.t)
        end
    end
end

@testset "time_filter on stiff problems" begin
    function rober!(du, u, p, t)
        du[1] = -0.04 * u[1] + 1.0e4 * u[2] * u[3]
        du[2] = 0.04 * u[1] - 1.0e4 * u[2] * u[3] - 3.0e7 * u[2]^2
        du[3] = 3.0e7 * u[2]^2
        return nothing
    end
    rober = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))

    reference = solve(rober, FBDF(), abstol = 1.0e-12, reltol = 1.0e-12)
    sol = solve(rober, FBDF(time_filter = true), abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test isapprox(sol.u[end], reference.u[end], rtol = 1.0e-5)

    # mass matrix DAE form: the filter's residual-based error estimate does not
    # apply, so this exercises the fallback estimate
    function rober_mm!(du, u, p, t)
        du[1] = -0.04 * u[1] + 1.0e4 * u[2] * u[3]
        du[2] = 0.04 * u[1] - 3.0e7 * u[2]^2 - 1.0e4 * u[2] * u[3]
        du[3] = u[1] + u[2] + u[3] - 1.0
        return nothing
    end
    rober_dae = ODEProblem(
        ODEFunction(rober_mm!, mass_matrix = Diagonal([1.0, 1.0, 0.0])),
        [1.0, 0.0, 0.0], (0.0, 1.0e5)
    )
    sol = solve(rober_dae, FBDF(time_filter = true), abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test isapprox(sol.u[end], reference.u[end], rtol = 1.0e-5)
end

@testset "time_filter respects max_order" begin
    for mo in (2, 3, 4, 5)
        sol = solve(
            IIP_DECAY, FBDF(max_order = Val(mo), time_filter = true),
            abstol = 1.0e-6, reltol = 1.0e-6
        )
        @test SciMLBase.successful_retcode(sol)
        @test decay_error(sol, 1.0e-6) < 100
    end
end

@testset "time_filter defaults to off" begin
    @test FBDF().time_filter == false
    for prob in (IIP_DECAY, OOP_DECAY)
        plain = solve(prob, FBDF(), abstol = 1.0e-8, reltol = 1.0e-8)
        off = solve(prob, FBDF(time_filter = false), abstol = 1.0e-8, reltol = 1.0e-8)
        @test SciMLBase.successful_retcode(plain)
        @test plain.t == off.t
        @test plain.u == off.u
    end
end
