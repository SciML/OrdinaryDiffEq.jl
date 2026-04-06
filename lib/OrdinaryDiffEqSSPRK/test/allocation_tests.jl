using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqCore
using Test

@inline function step_void!(integrator)
    step!(integrator)
    return nothing
end

@testset "SSPRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 100.0))

    ssprk_solvers = [
        SSPRK53_2N2(), SSPRK22(), SSPRK53(), SSPRK63(), SSPRK83(), SSPRK43(), SSPRK432(),
        SSPRKMSVS32(),
        SSPRK54(), SSPRK53_2N1(), SSPRK104(), SSPRK932(), SSPRKMSVS43(), SSPRK73(),
        SSPRK53_H(),
        SSPRK33(), KYKSSPRK42(), KYK2014DGSSPRK_3S2(),
        pRRK22(kappa = 0.0), pRRK22(kappa = 1.0),
        pRRK33(kappa = 0.0), pRRK33(kappa = 1.0),
        pRRK54(kappa = 0.0), pRRK54(kappa = 1.0),
    ]

    @testset "$(typeof(solver)) runtime allocation check" for solver in ssprk_solvers
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

            @testset "$(typeof(solver)) AllocCheck static analysis" for solver in ssprk_solvers
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
