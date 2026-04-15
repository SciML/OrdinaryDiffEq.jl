using OrdinaryDiffEqStabilizedRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

@testset "StabilizedRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    # All stabilized RK methods have allocations in perform_step!
    all_solvers = [ROCK2(), ROCK4(), RKC(), ESERK4(), ESERK5(), SERK2()]

    @testset "StabilizedRK perform_step! Static Analysis" begin
        for solver in all_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                step!(integrator)

                cache = integrator.cache
                # ROCK2/4 and ESERK use very large SVector fields (up to SVector{50})
                # which cause AllocCheck's LLVM analysis to throw MethodError.
                # Wrap in try-catch to gracefully skip those.
                allocs = try
                    check_allocs(
                        OrdinaryDiffEqCore.perform_step!,
                        (typeof(integrator), typeof(cache))
                    )
                catch e
                    println(
                        "check_allocs threw for $(typeof(solver)): $(typeof(e)) (skipping)"
                    )
                    nothing
                end

                if allocs !== nothing
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
end
