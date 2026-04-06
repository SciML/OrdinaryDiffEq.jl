using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqCore
using Test

@inline function step_void!(integrator)
    step!(integrator)
    return nothing
end

@testset "Rosenbrock Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 100.0))

    rosenbrock_solvers = [
        Rosenbrock23(), Rosenbrock32(), RosShamp4(), Veldd4(), Velds4(), GRK4T(), GRK4A(),
        Rodas3(), Rodas23W(), Rodas3P(), Rodas4(), Rodas42(), Rodas4P(), Rodas4P2(),
        Rodas5(),
        Rodas5P(), Rodas5Pe(), Rodas5Pr(), Rodas6P(),
    ]

    @testset "$(typeof(solver)) runtime allocation check" for solver in rosenbrock_solvers
        integrator = init(
            prob, solver, dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )

        for _ in 1:50
            step_void!(integrator)
        end

        GC.gc()
        for _ in 1:10
            step_void!(integrator)
        end

        GC.gc()
        allocs = 0
        for _ in 1:10
            allocs += @allocated step_void!(integrator)
        end
        @test allocs == 0
    end

    if !isempty(VERSION.prerelease)
        @info "Skipping AllocCheck static analysis on prerelease Julia"
    else
        try
            @eval using AllocCheck

            @testset "$(typeof(solver)) AllocCheck static analysis" for solver in rosenbrock_solvers
                fs_prob = ODEProblem{true, SciMLBase.FullSpecialize}(
                    simple_system!, [1.0, 1.0], (0.0, 1.0)
                )
                integrator = init(
                    fs_prob, solver, dt = 0.1, save_everystep = false,
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
                        "✓ $(typeof(solver)) perform_step! appears allocation-free with AllocCheck"
                    )
                end
            end
        catch e
            @info "AllocCheck not available, skipping static analysis" exception = e
        end
    end
end
