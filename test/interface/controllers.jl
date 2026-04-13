using OrdinaryDiffEq, OrdinaryDiffEqCore, Test

f(u, p, t) = zero(u)
prob = ODEProblem(f, 1.0, (0.0, 10.0))

# Shouldn't error
sol = solve(prob, Tsit5(), controller = IController())
@test sol.retcode == ReturnCode.Success

sol = solve(prob, Tsit5(), controller = PIController(7 // 50, 2 // 25))
@test sol.retcode == ReturnCode.Success

sol = solve(prob, Tsit5(), controller = PIDController(0.7, -0.4))
@test sol.retcode == ReturnCode.Success

# OrdinaryDiffEq.jl#1703
# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1703
prob = ODEProblem((du, u, p, t) -> du[1] = 1, [0.0], (0.0, 5.0))
sol = solve(prob, RDPK3SpFSAL49())
@test sol.retcode == ReturnCode.Success

# DifferentialEquations.jl#299
# Test that qmax_first_step allows larger step size growth on the first step.
# The initial dt from the automatic step size algorithm is only approximate,
# so the solver should allow larger step size ratios on the first step
# (matching Sundials CVODE behavior: 10^4 on first step vs 10 on later steps).
@testset "qmax_first_step" begin
    # Exponential decay - the auto dt init may pick a very conservative initial dt
    prob_exp = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 100.0))

    # Custom qmax_first_step via New* controller API
    ctrl = OrdinaryDiffEqCore.NewPIController(Tsit5(), qmax_first_step = 500)
    @test ctrl.qmax_first_step == 500.0

    # Verify that qmax_first_step=10000 still produces correct solutions
    for alg in (Tsit5(), Vern7(), RDPK3SpFSAL49())
        sol_default = solve(prob_exp, alg)
        @test sol_default.retcode == ReturnCode.Success
        @test isapprox(sol_default[end], exp(-100.0), atol = 1.0e-6)

        # With a very restrictive qmax_first_step (same as normal qmax),
        # the solution should still be correct
        ctrl_restricted = OrdinaryDiffEqCore.NewPIController(alg, qmax_first_step = 10)
        sol_restricted = solve(prob_exp, alg, controller = ctrl_restricted)
        @test sol_restricted.retcode == ReturnCode.Success
        @test isapprox(sol_restricted[end], exp(-100.0), atol = 1.0e-6)
    end

    # Test get_current_qmax: default 10000 on first step for all controllers
    integrator = init(prob_exp, Tsit5())
    @test integrator.success_iter == 0
    @test OrdinaryDiffEqCore.get_current_qmax(integrator, 10.0) == 10000.0
    # After a successful step, normal qmax is used
    integrator.success_iter = 1
    @test OrdinaryDiffEqCore.get_current_qmax(integrator, 10.0) == 10.0

    # Test get_current_qmax with New* controller (custom qmax_first_step)
    ctrl_new = OrdinaryDiffEqCore.NewPIController(Tsit5(), qmax_first_step = 500)
    integrator_new = init(prob_exp, Tsit5(), controller = ctrl_new)
    @test integrator_new.success_iter == 0
    @test OrdinaryDiffEqCore.get_current_qmax(integrator_new, 10.0) == 500.0
    integrator_new.success_iter = 1
    @test OrdinaryDiffEqCore.get_current_qmax(integrator_new, 10.0) == 10.0

    # Test that all New* controllers have qmax_first_step field
    ctrl_i = OrdinaryDiffEqCore.NewIController(Tsit5(), qmax_first_step = 200)
    @test ctrl_i.qmax_first_step == 200.0

    ctrl_pred = OrdinaryDiffEqCore.NewPredictiveController(Tsit5(), qmax_first_step = 300)
    @test ctrl_pred.qmax_first_step == 300.0
end
