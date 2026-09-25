using OrdinaryDiffEqSDC
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using Test

# Runtime measurement rather than AllocCheck: the static pass reports the
# `LazyString` sites behind `@SciMLMessage` in `OrdinaryDiffEqNonlinearSolve`,
# which a converging step never reaches.
@testset "SDC Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
        return nothing
    end

    # FullSpecialize avoids FunctionWrappers dynamic dispatch noise.
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    sdc_solvers = [
        SDC(),
        SDC(num_nodes = 4, num_sweeps = 4, sweeper = SDCSweeper.LU),
        SDC(num_nodes = 3, num_sweeps = 2, sweeper = SDCSweeper.FE),
        SDC(
            num_nodes = 3, num_sweeps = 3, quad_type = SDCQuadrature.Lobatto,
            step_update = SDCStepUpdate.LastNode
        ),
        SDC(num_nodes = 3, num_sweeps = 3, sweeper = SDCSweeper.BEpar),
    ]

    for solver in sdc_solvers
        @testset "$(nameof(typeof(solver))) $(solver.sweeper) M=$(solver.num_nodes) perform_step!" begin
            integrator = init(
                prob, solver, dt = 0.1, save_everystep = false, adaptive = false
            )
            cache = integrator.cache
            OrdinaryDiffEqCore.perform_step!(integrator, cache)
            @test @allocated(OrdinaryDiffEqCore.perform_step!(integrator, cache)) == 0
        end
    end
end
