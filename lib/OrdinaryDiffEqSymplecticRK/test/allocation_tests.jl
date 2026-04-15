using OrdinaryDiffEqSymplecticRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, DynamicalODEFunction
using AllocCheck
using Test

@testset "SymplecticRK Allocation Tests" begin
    # Simple harmonic oscillator: dv = -u, du = v
    function f1!(dv, v, u, p, t)
        dv .= -u
    end
    function f2!(du, v, u, p, t)
        du .= v
    end

    v0 = [1.0, 0.0]
    u0 = [0.0, 1.0]

    # FullSpecialize on DynamicalODEFunction avoids FunctionWrappers noise
    ff = DynamicalODEFunction{true, FullSpecialize}(f1!, f2!)
    prob = DynamicalODEProblem(ff, v0, u0, (0.0, 1.0))

    symplectic_solvers = [
        SymplecticEuler(), VelocityVerlet(), VerletLeapfrog(), LeapfrogDriftKickDrift(),
        PseudoVerletLeapfrog(), McAte2(), Ruth3(), McAte3(), CandyRoz4(),
        McAte4(), McAte42(), McAte5(), CalvoSanz4(), Yoshida6(),
        KahanLi6(), McAte8(), KahanLi8(), SofSpa10(),
    ]

    @testset "SymplecticRK perform_step! Static Analysis" begin
        for solver in symplectic_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.01, save_everystep = false, adaptive = false
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
