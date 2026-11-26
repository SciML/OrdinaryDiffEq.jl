using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqRosenbrock, Test

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/4801
# u1(t) = 1/(1 - t) exactly; u2 is a clock (u2 == t for any consistent method).
f_blowup!(du, u, p, t) = (du[1] = u[1]^2; du[2] = 1.0; nothing)
f_blowup(u, p, t) = [u[1]^2, 1.0]

function f_rober!(du, u, p, t)
    y1, y2, y3 = u
    du[1] = -0.04y1 + 1.0e4 * y2 * y3
    du[2] = 0.04y1 - 1.0e4 * y2 * y3 - 3.0e7 * y2^2
    du[3] = 3.0e7 * y2^2
    return nothing
end
f_rober(u, p, t) = (du = similar(u); f_rober!(du, u, p, t); du)

blowup_probs = (
    ODEProblem(f_blowup!, [1.0, 0.0], (0.0, 0.9)),
    ODEProblem(f_blowup, [1.0, 0.0], (0.0, 0.9)),
)
rober_probs = (
    ODEProblem(f_rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5)),
    ODEProblem(f_rober, [1.0, 0.0, 0.0], (0.0, 1.0e5)),
)
rober_ref = solve(rober_probs[1], Rodas5P(); reltol = 1.0e-12, abstol = 1.0e-14).u[end]

pid_controllers(alg) = (
    OrdinaryDiffEqCore.PIDController(0.7, -0.4),
    OrdinaryDiffEqCore.PIDController(alg),
)

# A run is stalled if it needs far more steps than the default controller (master took
# ~1e6 accepts at dt ~ 1e-8 where the default controller takes ~1e3).
const MAXITERS = 100_000

@testset "PIDController with BDF-family methods" begin
    @testset "PIDController acceptance agrees with the step commit" begin
        for alg in (QNDF1(), QNDF2(), QNDF(), FBDF(), ABDF2()), prob in blowup_probs,
                controller in pid_controllers(alg)
            integ = init(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6, controller)
            accepted_with_large_EEst = 0
            while integ.t < prob.tspan[2] && integ.iter < MAXITERS
                step!(integ)
                integ.accept_step &&
                    OrdinaryDiffEqCore.get_EEst(integ) > 1 &&
                    (accepted_with_large_EEst += 1)
            end
            @test accepted_with_large_EEst == 0
            @test integ.t == prob.tspan[2]
        end
    end

    @testset "$(nameof(typeof(alg))) + PIDController on the blowup problem" for alg in (
            QNDF1(), QNDF2(), QNDF(), FBDF(), ABDF2(),
        )
        for prob in blowup_probs
            default = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6)
            for controller in pid_controllers(alg)
                sol = solve(
                    prob, alg; reltol = 1.0e-6, abstol = 1.0e-6, controller, maxiters = MAXITERS
                )
                @test sol.retcode == ReturnCode.Success
                @test sol.t[end] == 0.9
                # The clock is exact up to the Newton tolerance, so a lost step shows here.
                @test sol.u[end][2] ≈ 0.9 atol = 1.0e-6
                # Exact value 1/(1 - 0.9) = 10; the default controller is within 1.1% here.
                @test sol.u[end][1] ≈ 10 rtol = 0.02
                @test sol.stats.naccept <= 2 * default.stats.naccept
            end
        end
    end

    @testset "$(nameof(typeof(alg))) + PIDController on Robertson" for alg in (
            QNDF1(), QNDF2(), QNDF(), FBDF(),
        )
        for prob in rober_probs
            default = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-10)
            for controller in pid_controllers(alg)
                sol = solve(
                    prob, alg; reltol = 1.0e-6, abstol = 1.0e-10, controller, maxiters = MAXITERS
                )
                @test sol.retcode == ReturnCode.Success
                @test sol.t[end] == 1.0e5
                # The default controller reaches y1, y3 within 9e-4 relative of the Rodas5P
                # reference with QNDF1, the least accurate method here.
                @test sol.u[end][1] ≈ rober_ref[1] rtol = 5.0e-3
                @test sol.u[end][3] ≈ rober_ref[3] rtol = 5.0e-3
                @test sum(sol.u[end]) ≈ 1 atol = 1.0e-8
                @test sol.stats.naccept <= 2 * default.stats.naccept
            end
        end
    end
end
