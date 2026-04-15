using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize
using AllocCheck
using Test

@testset "LowStorageRK Allocation Tests" begin
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end

    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise
    prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))

    # 2N low-storage methods
    # RDPK3 family excluded: PID controller init bug (MVector{3} from Bool flags)
    lsrk_2n_solvers = [
        ORK256(), CarpenterKennedy2N54(), SHLDDRK64(), HSLDDRK64(),
        DGLDDRK73_C(), DGLDDRK84_C(), DGLDDRK84_F(),
        NDBLSRK124(), NDBLSRK134(), NDBLSRK144(),
        CFRLDDRK64(), TSLDDRK74(),
        CKLLSRK43_2(), CKLLSRK54_3C(), CKLLSRK95_4S(), CKLLSRK95_4C(), CKLLSRK95_4M(),
        CKLLSRK54_3C_3R(), CKLLSRK54_3M_3R(), CKLLSRK54_3N_3R(),
        CKLLSRK85_4C_3R(), CKLLSRK85_4M_3R(), CKLLSRK85_4P_3R(),
        CKLLSRK54_3N_4R(), CKLLSRK54_3M_4R(), CKLLSRK65_4M_4R(),
        CKLLSRK85_4FM_4R(), CKLLSRK75_4M_5R(),
        ParsaniKetchesonDeconinck3S32(), ParsaniKetchesonDeconinck3S82(),
        ParsaniKetchesonDeconinck3S53(), ParsaniKetchesonDeconinck3S173(),
        ParsaniKetchesonDeconinck3S94(), ParsaniKetchesonDeconinck3S184(),
        ParsaniKetchesonDeconinck3S105(), ParsaniKetchesonDeconinck3S205(),
        RK46NL(), SHLDDRK_2N(), SHLDDRK52(),
    ]

    @testset "LowStorageRK perform_step! Static Analysis" begin
        for solver in lsrk_2n_solvers
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
    end
end
