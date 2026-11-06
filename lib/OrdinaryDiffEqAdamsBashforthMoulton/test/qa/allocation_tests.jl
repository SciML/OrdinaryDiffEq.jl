using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

@testset "AdamsBashforthMoulton Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 10.0))

    # Solvers confirmed allocation-free in `perform_step!`. Any regression
    # here fails the test loud.
    # AllocCheck reports broadcast-alias sites for ABM32 outside Julia 1.12 bounds=auto.
    abm32_alloc_free = VERSION >= v"1.12" && Base.JLOptions().check_bounds == 0
    ab_solvers_alloc_free = abm32_alloc_free ? [AB3(), ABM32()] : [AB3()]

    # Solvers with known allocation sites in `perform_step!`. Marked
    # `broken = true` so the suite stays green while tracking the
    # regression. When a solver is fixed, AllocCheck will report zero
    # sites, the test flips to "Unexpected Pass", and the entry should
    # move up to the `*_alloc_free` list above.
    ab_solvers_known_broken = abm32_alloc_free ?
        [AB4(), AB5(), ABM43(), ABM54()] :
        [AB4(), AB5(), ABM32(), ABM43(), ABM54()]
    vcab_solvers_known_broken = [
        VCAB3(), VCAB4(), VCAB5(),
        VCABM3(), VCABM4(), VCABM5(), VCABM(),
    ]

    @testset "ABM perform_step! Static Analysis" begin
        # Fixed-step AB / ABM — alloc-free enforced
        for solver in ab_solvers_alloc_free
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false, adaptive = false
                )
                # Multistep methods need history: advance past startup steps
                for _ in 1:5
                    step!(integrator)
                end

                cache = integrator.cache
                allocs = check_allocs(
                    OrdinaryDiffEqCore.perform_step!,
                    (typeof(integrator), typeof(cache))
                )

                @test length(allocs) == 0

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

        # Fixed-step AB / ABM — known allocation regressions
        for solver in ab_solvers_known_broken
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false, adaptive = false
                )
                for _ in 1:5
                    step!(integrator)
                end

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

        # Variable-coefficient Adams — all currently have allocation sites
        for solver in vcab_solvers_known_broken
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                for _ in 1:5
                    step!(integrator)
                end

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

    @testset "step! runtime allocations" begin
        # AllocCheck over-reports broadcast alias copies, so the steady
        # step is also required to be allocation-free at runtime. step!
        # itself carries a small fixed allocation on Julia LTS that is
        # absent on 1.12, so each solver is compared against the AB3
        # baseline instead of zero.
        ref_integrator = init(
            prob, AB3(), dt = 0.1, save_everystep = false, adaptive = false
        )
        for _ in 1:10
            step!(ref_integrator)
        end
        ref_bytes = @allocated step!(ref_integrator)

        for solver in [
                ABM32(), ABM43(), ABM54(),
                VCAB3(), VCAB4(), VCAB5(),
                VCABM3(), VCABM4(), VCABM5(), VCABM(),
            ]
            @testset "$(typeof(solver)) step! is allocation-free" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
                    adaptive = OrdinaryDiffEqCore.isadaptive(solver),
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                # Multistep methods need history: advance past startup steps
                for _ in 1:10
                    step!(integrator)
                end
                @test @allocated(step!(integrator)) == ref_bytes
            end
        end
    end
end
