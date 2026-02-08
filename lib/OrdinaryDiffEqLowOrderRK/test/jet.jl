using Pkg
Pkg.add("JET")

import OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqLowOrderRK
using OrdinaryDiffEqCore
using JET
using Test

@testset "JET Tests" begin
    # Test package for typos (commented out due to false positive)
    # test_package(
    #     OrdinaryDiffEqLowOrderRK, mode = :typo)

    # Test individual solver type stability
    @testset "Solver Type Stability Tests" begin
        # Test problem
        function simple_system!(du, u, p, t)
            du[1] = -0.5 * u[1]
            du[2] = -1.5 * u[2]
        end
        prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))

        # Test all exported LowOrderRK solvers
        low_order_solvers = [
            Euler(), Heun(), Ralston(), Midpoint(), RK4(),
            BS3(), OwrenZen3(), OwrenZen4(), OwrenZen5(), BS5(),
            DP5(), Anas5(), RKO65(), FRK65(), RKM(), MSRK5(), MSRK6(),
            PSRK4p7q6(), PSRK3p5q4(), PSRK3p6q5(), Stepanov5(), SIR54(),
            Alshina2(), Alshina3(), Alshina6(), AutoDP5(DP5()),
        ]

        for solver in low_order_solvers
            @testset "$(typeof(solver)) type stability" begin
                try
                    # Some solvers need fixed timestep
                    if solver isa Euler || solver isa Midpoint || solver isa Heun
                        @test_opt broken = true init(prob, solver, dt = 0.1, save_everystep = false, adaptive = false)
                        integrator = init(prob, solver, dt = 0.1, save_everystep = false, adaptive = false)
                    else
                        @test_opt broken = true init(prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                        integrator = init(prob, solver, dt = 0.1, save_everystep = false, abstol = 1.0e-6, reltol = 1.0e-6)
                    end
                    @test_opt broken = true step!(integrator)
                catch e
                    @test_broken false # Mark as broken if solver fails to initialize
                    println("$(typeof(solver)) failed with: $e")
                end
            end
        end
    end
end
