using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

@testset "Rosenbrock Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    rosenbrock_solvers = [
        Rosenbrock23(), Rosenbrock32(), RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
        Rodas3(), Rodas23W(), Rodas3P(), Rodas4(), Rodas42(), Rodas4P(), Rodas4P2(), Rodas5(),
        Rodas5P(), Rodas5Pe(), Rodas5Pr(), Rodas6P(),
    ]

    @testset "Rosenbrock perform_step! Static Analysis" begin
        for solver in rosenbrock_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                step!(integrator)

                cache = integrator.cache
                allocs = check_allocs(
                    OrdinaryDiffEqCore.perform_step!,
                    (typeof(integrator), typeof(cache))
                )

                @test length(allocs) == 0 broken = true

                if length(allocs) > 0
                    println(
                        "AllocCheck found $(length(allocs)) allocation sites in $(typeof(solver)) perform_step!"
                    )
                else
                    println(
                        "$(typeof(solver)) perform_step! appears allocation-free with AllocCheck"
                    )
                end
            end
        end
    end
end
