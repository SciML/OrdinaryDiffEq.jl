using OrdinaryDiffEqPRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

@testset "PRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    @testset "KuttaPRK2p5 perform_step! Static Analysis" begin
        integrator = init(
            prob, KuttaPRK2p5(), dt = 0.1, save_everystep = false,
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
                "AllocCheck found $(length(allocs)) allocation sites in KuttaPRK2p5 perform_step!"
            )
        else
            println("KuttaPRK2p5 perform_step! appears allocation-free with AllocCheck")
        end
    end
end
