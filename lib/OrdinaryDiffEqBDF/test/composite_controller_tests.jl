using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, Test
using OrdinaryDiffEqCore: CompositeAlgorithm, accept_step_controller, get_EEst, set_EEst!

@testset "BDF composite error acceptance" begin
    prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    for alg in (FBDF(), QNDF())
        composite = CompositeAlgorithm((alg,), _ -> 1)
        integrator = init(prob, composite)
        for estimate in (0.5, 1.0, 2.0, NaN)
            set_EEst!(integrator, estimate)
            @test accept_step_controller(integrator, composite) == (estimate <= 1)
        end
    end
end

@testset "FBDF composite step rejection" begin
    rates = [1.0, 10.0, 100.0]
    f(u, p, t) = -p .* u
    f!(du, u, p, t) = (du .= -p .* u; nothing)
    for rhs in (f, f!)
        prob = ODEProblem(rhs, ones(3), (0.0, 1.0), rates)
        reference = init(prob, FBDF(); dt = 0.1, abstol = 1.0e-12, reltol = 1.0e-9)
        solve!(reference)
        @test reference.sol.retcode == ReturnCode.Success
        @test reference.stats.nreject > 0
        for alg in (
                CompositeAlgorithm((FBDF(),), _ -> 1),
                CompositeAlgorithm((QNDF(), FBDF()), _ -> 2),
            )
            integrator = init(prob, alg; dt = 0.1, abstol = 1.0e-12, reltol = 1.0e-9)
            accepted_over_1 = 0
            for _ in integrator
                accepted_over_1 += get_EEst(integrator) > 1
            end
            @test integrator.sol.retcode == ReturnCode.Success
            @test accepted_over_1 == 0
            @test integrator.stats.nreject == reference.stats.nreject
            @test integrator.sol.t == reference.sol.t
            @test integrator.sol.u == reference.sol.u
        end
    end
end
