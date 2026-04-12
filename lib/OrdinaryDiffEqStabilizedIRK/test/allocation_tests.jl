using OrdinaryDiffEqStabilizedIRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, SplitFunction, ODEFunction
using AllocCheck
using LinearAlgebra
using Test

@testset "StabilizedIRK Allocation Tests" begin
    A = [-10.0 0.5; 0.0 -2.0]
    function f1!(du, u, p, t)
        mul!(du, A, u)
    end
    function f2!(du, u, p, t)
        du[1] = 0.0
        du[2] = 0.0
    end

    ff = SplitFunction{true}(
        ODEFunction{true, FullSpecialize}(f1!),
        ODEFunction{true, FullSpecialize}(f2!)
    )
    prob = SplitODEProblem(ff, [1.0, 1.0], (0.0, 1.0))

    @testset "IRKC perform_step! Static Analysis" begin
        integrator = init(
            prob, IRKC(), dt = 0.1, save_everystep = false,
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
                "AllocCheck found $(length(allocs)) allocation sites in IRKC perform_step!"
            )
        else
            println("IRKC perform_step! appears allocation-free with AllocCheck")
        end
    end
end
