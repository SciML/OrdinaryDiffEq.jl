using OrdinaryDiffEqBDF
using SciMLBase
using LinearAlgebra
using Test

exp_decay!(du, u, p, t) = (du[1] = -u[1]; nothing)
exp_decay(u, p, t) = -u

const IIP_DECAY = ODEProblem(exp_decay!, [1.0], (0.0, 10.0))
const OOP_DECAY = ODEProblem(exp_decay, 1.0, (0.0, 10.0))

function decay_error(sol, tol)
    return maximum(
        abs(only(u) - exp(-t)) / (tol + tol * exp(-t)) for (u, t) in zip(sol.u, sol.t)
    )
end

@testset "NordsieckBDF time_filter defaults to off" begin
    @test NordsieckBDF().time_filter == false
    @test NordsieckBDF(time_filter = true).time_filter == true
end

@testset "NordsieckBDF time_filter solves accurately" begin
    for (name, prob) in (("in-place", IIP_DECAY), ("out-of-place", OOP_DECAY))
        @testset "$name" begin
            for tol in (1.0e-6, 1.0e-8, 1.0e-10)
                sol = solve(
                    prob, NordsieckBDF(time_filter = true), abstol = tol, reltol = tol
                )
                @test SciMLBase.successful_retcode(sol)
                @test decay_error(sol, tol) < 100
            end
        end
    end
end

@testset "NordsieckBDF time_filter does not need more steps on a smooth problem" begin
    for prob in (IIP_DECAY, OOP_DECAY)
        for tol in (1.0e-6, 1.0e-8)
            filtered = solve(
                prob, NordsieckBDF(time_filter = true), abstol = tol, reltol = tol
            )
            plain = solve(prob, NordsieckBDF(), abstol = tol, reltol = tol)
            @test SciMLBase.successful_retcode(filtered)
            @test length(filtered.t) <= 1.15 * length(plain.t)
        end
    end
end

@testset "NordsieckBDF time_filter on stiff problems" begin
    function rober!(du, u, p, t)
        du[1] = -0.04 * u[1] + 1.0e4 * u[2] * u[3]
        du[2] = 0.04 * u[1] - 1.0e4 * u[2] * u[3] - 3.0e7 * u[2]^2
        du[3] = 3.0e7 * u[2]^2
        return nothing
    end
    rober = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5))

    reference = solve(rober, NordsieckBDF(), abstol = 1.0e-12, reltol = 1.0e-12)
    sol = solve(rober, NordsieckBDF(time_filter = true), abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
    @test isapprox(sol.u[end], reference.u[end], rtol = 1.0e-5)

    for iip in (false, true)
        dae_f(u, p, t) = [-u[1], u[2] - u[1]^2]
        dae_f!(du, u, p, t) = (du .= dae_f(u, p, t); nothing)
        fun = iip ? ODEFunction(dae_f!; mass_matrix = Diagonal([1.0, 0.0])) :
            ODEFunction(dae_f; mass_matrix = Diagonal([1.0, 0.0]))
        dae = ODEProblem(fun, [1.0, 1.0], (0.0, 1.0))
        @test_throws ArgumentError init(dae, NordsieckBDF(time_filter = true))
        @test SciMLBase.successful_retcode(solve(dae, NordsieckBDF(time_filter = false)))
    end
end

@testset "NordsieckBDF time_filter respects max_order" begin
    for mo in (2, 3, 4, 5)
        sol = solve(
            IIP_DECAY, NordsieckBDF(max_order = Val(mo), time_filter = true),
            abstol = 1.0e-6, reltol = 1.0e-6
        )
        @test SciMLBase.successful_retcode(sol)
        @test decay_error(sol, 1.0e-6) < 100
    end
end

@testset "NordsieckBDF time_filter fixed-step oscillator" begin
    function osc!(du, u, p, t)
        ω = p
        du[1] = -ω * u[2]
        du[2] = ω * u[1]
        return nothing
    end
    function osc_l2(sol, ω)
        err2 = 0.0
        for (u, t) in zip(sol.u, sol.t)
            err2 += sum(abs2, u .- (cos(ω * t), sin(ω * t)))
        end
        return sqrt(err2 / length(sol.t))
    end
    # Mild oscillator: filter must not spoil an already-accurate solve.
    osc_mild = ODEProblem(osc!, [1.0, 0.0], (0.0, 2.0), 10.0)
    plain = solve(osc_mild, NordsieckBDF(); adaptive = false, dt = 1.0e-3, dense = false)
    filt = solve(
        osc_mild, NordsieckBDF(time_filter = true);
        adaptive = false, dt = 1.0e-3, dense = false
    )
    @test SciMLBase.successful_retcode(plain)
    @test SciMLBase.successful_retcode(filt)
    @test osc_l2(filt, 10.0) < 2 * osc_l2(plain, 10.0)

    # Stiffer imaginary eigenvalue / larger dt: filter should improve accuracy.
    osc_hard = ODEProblem(osc!, [1.0, 0.0], (0.0, 2.0), 100.0)
    plain_h = solve(osc_hard, NordsieckBDF(); adaptive = false, dt = 1.0e-3, dense = false)
    filt_h = solve(
        osc_hard, NordsieckBDF(time_filter = true);
        adaptive = false, dt = 1.0e-3, dense = false
    )
    @test SciMLBase.successful_retcode(plain_h)
    @test SciMLBase.successful_retcode(filt_h)
    @test osc_l2(filt_h, 100.0) <= osc_l2(plain_h, 100.0)
end
