using OrdinaryDiffEqLinear
using OrdinaryDiffEqCore
using AllocCheck
using LinearAlgebra
using Test

@testset "Linear Allocation Tests" begin
    A = MatrixOperator([-1.0 0.5; 0.0 -2.0])
    prob = ODEProblem(A, [1.0, 1.0], (0.0, 1.0))

    # CayleyEuler excluded: requires matrix-valued state + SplitODEProblem
    linear_solvers = [
        MagnusMidpoint(), LieEuler(), RKMK2(), RKMK4(), LieRK4(),
        CG2(), CG3(), CG4a(), LinearExponential(krylov = :off), MagnusAdapt4(),
    ]

    @testset "Linear perform_step! Static Analysis" begin
        for solver in linear_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false, adaptive = false
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
