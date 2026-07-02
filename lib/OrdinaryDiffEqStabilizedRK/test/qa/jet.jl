import OrdinaryDiffEqStabilizedRK
using OrdinaryDiffEqStabilizedRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using JET
using Test

@testset "JET Tests" begin
    test_package(
        OrdinaryDiffEqStabilizedRK, target_modules = (OrdinaryDiffEqStabilizedRK,), mode = :typo
    )

    @testset "perform_step! type stability" begin
        function stiff_system!(du, u, p, t)
            du[1] = -20.0 * u[1]
            du[2] = -10.0 * u[2]
        end
        prob = ODEProblem{true, FullSpecialize}(stiff_system!, [1.0, 1.0], (0.0, 1.0))
        oop_prob = ODEProblem{false, FullSpecialize}((u, p, t) -> -20.0 * u, 1.0, (0.0, 1.0))

        all_solvers = [
            ROCK2(), ROCK4(), RKC(), RKMC2(), ESERK4(), ESERK5(), SERK2(),
            TSRKC2(), TSRKC3(), RKL1(), RKL2(), RKG1(), RKG2(),
        ]

        for solver in all_solvers
            @testset "$(nameof(typeof(solver))) in-place" begin
                integrator = init(
                    prob, solver, dt = 0.01, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                step!(integrator)
                @test_opt target_modules = (OrdinaryDiffEqStabilizedRK,) OrdinaryDiffEqCore.perform_step!(
                    integrator, integrator.cache
                )
            end
            @testset "$(nameof(typeof(solver))) out-of-place" begin
                integrator = init(
                    oop_prob, solver, dt = 0.01, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                step!(integrator)
                @test_opt target_modules = (OrdinaryDiffEqStabilizedRK,) OrdinaryDiffEqCore.perform_step!(
                    integrator, integrator.cache
                )
            end
        end
    end
end
