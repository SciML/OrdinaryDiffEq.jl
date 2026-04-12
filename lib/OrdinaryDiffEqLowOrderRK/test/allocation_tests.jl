using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

"""
Allocation tests for OrdinaryDiffEqLowOrderRK solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.
Uses FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.
"""

@testset "LowOrderRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    low_order_solvers = [
        Euler(), Heun(), Ralston(), Midpoint(), RK4(),
        BS3(), OwrenZen3(), OwrenZen4(), OwrenZen5(), BS5(),
        DP5(), Anas5(), RKO65(), FRK65(), RKM(), MSRK5(), MSRK6(),
        PSRK4p7q6(), PSRK3p5q4(), PSRK3p6q5(), SIR54(),
        Alshina2(), Alshina3(), Alshina6(), AutoDP5(DP5()),
    ]

    # Stepanov5 has real allocations in perform_step!
    low_order_broken_solvers = [Stepanov5()]

    @testset "LowOrderRK perform_step! Static Analysis" begin
        for solver in low_order_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                if solver isa Euler || solver isa Midpoint || solver isa Heun
                    integrator = init(
                        prob, solver, dt = 0.1, save_everystep = false, adaptive = false
                    )
                else
                    integrator = init(
                        prob, solver, dt = 0.1, save_everystep = false,
                        abstol = 1.0e-6, reltol = 1.0e-6
                    )
                end
                step!(integrator)

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

        for solver in low_order_broken_solvers
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
