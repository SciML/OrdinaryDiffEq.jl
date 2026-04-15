using OrdinaryDiffEqSIMDRK
using OrdinaryDiffEqCore
using AllocCheck
using Test

@testset "SIMDRRK Allocation Tests" begin
    # OOP (out-of-place) problem — SIMDRRK only supports non-IIP
    function simple_system(u, p, t)
        return [-0.5 * u[1], -1.5 * u[2]]
    end

    prob = ODEProblem(simple_system, [1.0, 1.0], (0.0, 1.0))

    # RK6v4 excluded: requires SIMD packed-vector types, incompatible with Vector{Float64}
    simd_solvers = [MER5v2(), MER6v2()]

    @testset "SIMDRRK perform_step! Static Analysis" begin
        for solver in simd_solvers
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
