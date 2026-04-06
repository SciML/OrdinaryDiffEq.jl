using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqCore
using Test

@inline function step_void!(integrator)
    step!(integrator)
    return nothing
end

@testset "LowOrderRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 100.0))

    # Adaptive solvers
    adaptive_solvers = [
        BS3(), OwrenZen3(), OwrenZen4(), OwrenZen5(), BS5(),
        DP5(), Anas5(), RKO65(), FRK65(), RKM(), MSRK5(), MSRK6(),
        PSRK4p7q6(), PSRK3p5q4(), PSRK3p6q5(), SIR54(),
        Alshina2(), Alshina3(), Alshina6(),
    ]

    # Solvers with known allocations (tested with @test_broken)
    known_allocating_solvers = [Stepanov5()]

    # Fixed-timestep solvers
    fixed_dt_solvers = [Euler(), Heun(), Ralston(), Midpoint(), RK4()]

    # Auto-switching solvers (may allocate due to algorithm switching)
    auto_solvers = [AutoDP5(DP5())]

    @testset "$(typeof(solver)) runtime allocation check" for solver in adaptive_solvers
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

    @testset "$(typeof(solver)) runtime allocation check" for solver in known_allocating_solvers
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
        @test allocs == 0 broken = true
    end

    @testset "$(typeof(solver)) runtime allocation check" for solver in fixed_dt_solvers
        integrator = init(
            prob, solver, dt = 0.1, save_everystep = false,
            adaptive = false
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

    @testset "$(typeof(solver)) runtime allocation check" for solver in auto_solvers
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

            all_solvers = vcat(adaptive_solvers, fixed_dt_solvers)

            @testset "$(typeof(solver)) AllocCheck static analysis" for solver in all_solvers
                fs_prob = ODEProblem{true, SciMLBase.FullSpecialize}(
                    simple_system!, [1.0, 1.0], (0.0, 1.0)
                )
                if solver isa Euler || solver isa Midpoint || solver isa Heun ||
                        solver isa Ralston || solver isa RK4
                    integrator = init(
                        fs_prob, solver, dt = 0.1, save_everystep = false,
                        adaptive = false
                    )
                else
                    integrator = init(
                        fs_prob, solver, dt = 0.1, save_everystep = false,
                        abstol = 1.0e-6, reltol = 1.0e-6
                    )
                end
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
