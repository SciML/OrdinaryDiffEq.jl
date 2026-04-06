using OrdinaryDiffEqVerner
using OrdinaryDiffEqCore
using Test

@inline function step_void!(integrator)
    step!(integrator)
    return nothing
end

@testset "Verner Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 100.0))

    # Core Verner solvers (should be allocation-free)
    verner_solvers = [Vern6(), Vern7(), Vern8(), Vern9()]

    # AutoVern solvers involve algorithm switching which may allocate
    auto_verner_solvers = [
        AutoVern6(Vern6()), AutoVern7(Vern7()),
        AutoVern8(Vern8()), AutoVern9(Vern9()),
    ]

    @testset "$(typeof(solver)) runtime allocation check" for solver in verner_solvers
        integrator = init(
            prob, solver, dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )

        for _ in 1:10
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

    # AutoVern solvers - algorithm switching may allocate, use @test_broken
    @testset "$(typeof(solver)) runtime allocation check" for solver in auto_verner_solvers
        integrator = init(
            prob, solver, dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )

        for _ in 1:10
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

            @testset "$(typeof(solver)) AllocCheck static analysis" for solver in verner_solvers
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
