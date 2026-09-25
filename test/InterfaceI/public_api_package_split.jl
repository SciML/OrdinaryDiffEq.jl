using DiffEqCallbacks: FunctionCallingCallback
using OrdinaryDiffEqBDF: FBDF
using OrdinaryDiffEqLowOrderRK: DP5
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqSDIRK: TRBDF2
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqVerner: AutoVern9, Vern9
using RecursiveArrayTools: ArrayPartition
using SciMLBase: ODEProblem, ReturnCode, solve, successful_retcode, terminate!
using Test

@testset "Public API package splits" begin
    linear_prob = ODEProblem((u, p, t) -> p .* u, [1.0], (0.0, 1.0), -1.0)

    for (name, alg) in (
            ("Tsit5", Tsit5()),
            ("DP5", DP5()),
            ("Vern9", Vern9()),
            ("Rodas5P", Rodas5P()),
            ("TRBDF2", TRBDF2()),
            ("FBDF", FBDF()),
        )
        @testset "$name" begin
            sol = solve(linear_prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)

            @test successful_retcode(sol)
            @test sol.u[end][1] ≈ exp(-1) rtol = 1.0e-4
        end
    end

    function chart_exp_problem(u, p, t)
        return ArrayPartition(u.x[2], zero(u.x[2]))
    end

    function terminate_after_tstop(u, t, integrator)
        t >= 0.6 && terminate!(integrator)
        return u
    end

    u0 = ArrayPartition([0.0], [1.0])
    callback = FunctionCallingCallback(terminate_after_tstop; func_start = false)
    prob = ODEProblem(chart_exp_problem, u0, (0.0, 1.0); callback)

    @testset "Auto-switch callback and structured state" begin
        sol = solve(prob, AutoVern9(Rodas5P()); save_everystep = false, tstops = [0.6])

        @test sol.retcode === ReturnCode.Terminated
        @test sol.u[end] isa ArrayPartition
        @test sol.t[end] ≈ 0.6 atol = 1.0e-10
    end
end
