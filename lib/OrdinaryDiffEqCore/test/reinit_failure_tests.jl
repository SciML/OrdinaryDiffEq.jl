using OrdinaryDiffEqRosenbrock, Test

@testset "Reinit clears failed-step state" begin
    for (alg, gamma) in (
                (Rodas5P(), OrdinaryDiffEqRosenbrock.Rodas5PTableau(Float64, Float64).gamma),
                (Rosenbrock23(), OrdinaryDiffEqRosenbrock.Rosenbrock23Tableau(Float64).d),
                (Rosenbrock32(), OrdinaryDiffEqRosenbrock.Rosenbrock32Tableau(Float64).d),
            ), adaptive in (false, true), reset_dt in (false, true)
        dt = 0.125
        # W is exactly singular on the first step. Changing p below makes
        # the reinitialized problem a regular decay equation.
        prob = ODEProblem((u, p, t) -> p .* u, [1.0], (0.0, 0.5), inv(gamma * dt))
        kw = (; dt, adaptive, dtmin = dt, dtmax = dt, force_dtmin = true, maxiters = 25)
        integrator = init(prob, alg; kw...)
        @test solve!(integrator).retcode == ReturnCode.ConvergenceFailure
        @test integrator.last_stepfail

        integrator.p = -1.0
        reinit!(integrator, [2.0]; reset_dt)
        @test !integrator.force_stepfail
        @test !integrator.last_stepfail
        sol = solve!(integrator)
        reference = solve(remake(prob; p = -1.0, u0 = [2.0]), alg; kw...)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t[end] == prob.tspan[2]
        @test sol.u[end] ≈ reference.u[end]
    end
end

@testset "Reinit respects return-code options and callback initialization" begin
    alg = Rodas5P()
    gamma = OrdinaryDiffEqRosenbrock.Rodas5PTableau(Float64, Float64).gamma
    dt = 0.125
    prob = ODEProblem((u, p, t) -> p .* u, [1.0], (0.0, 0.5), inv(gamma * dt))
    for reinit_retcode in (false, true)
        integrator = init(prob, alg; dt, adaptive = false)
        @test solve!(integrator).retcode == ReturnCode.ConvergenceFailure
        integrator.p = -1.0
        reinit!(integrator, [2.0]; reinit_retcode)
        @test !integrator.last_stepfail
        @test integrator.sol.retcode == (reinit_retcode ? ReturnCode.Default : ReturnCode.ConvergenceFailure)
        sol = solve!(integrator)
        @test sol.retcode == (reinit_retcode ? ReturnCode.Success : ReturnCode.ConvergenceFailure)
    end

    initialize = (cb, u, t, integrator) -> (u[1] > 1.5 && terminate!(integrator))
    callback = DiscreteCallback((u, t, integrator) -> false, identity; initialize)
    integrator = init(prob, alg; dt, adaptive = false, callback)
    @test solve!(integrator).retcode == ReturnCode.ConvergenceFailure
    integrator.p = -1.0
    reinit!(integrator, [2.0])
    @test !integrator.last_stepfail
    @test integrator.sol.retcode == ReturnCode.Terminated
    @test solve!(integrator).retcode == ReturnCode.Terminated
end
