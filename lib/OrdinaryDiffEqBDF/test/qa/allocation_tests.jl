using OrdinaryDiffEqBDF
using OrdinaryDiffEqCore
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using SciMLBase: FullSpecialize, SplitFunction, ODEFunction, DAEFunction
using AllocCheck
using ADTypes: AutoForwardDiff, AutoSparse
using OrdinaryDiffEqDifferentiation: OrdinaryDiffEqDifferentiation
using SciMLBase: NoInit
using SparseArrays
using SparseConnectivityTracer: TracerSparsityDetector
using SparseMatrixColorings: GreedyColoringAlgorithm
using Test

@testset "BDF Allocation Tests" begin
    # Test problem for standard ODE BDF methods
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))

    # Split problem for IMEX methods
    function f1!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = 0.0
    end
    function f2!(du, u, p, t)
        du[1] = 0.0
        du[2] = -1.5 * u[2]
    end
    split_prob = SplitODEProblem(f1!, f2!, [1.0, 1.0], (0.0, 1.0))

    # DAE problem for DAE solvers
    function dae_f!(resid, du, u, p, t)
        resid[1] = -0.5 * u[1] + u[2] - du[1]
        resid[2] = u[1] - u[2] - du[2]
    end
    du0 = zeros(2)
    differential_vars = [true, false]
    dae_prob = DAEProblem(
        dae_f!, du0, [1.0, 1.0], (0.0, 1.0); differential_vars
    )

    # Runtime allocation tests run FIRST, before AllocCheck static analysis
    # which can invalidate compiled code and cause false positive allocations.

    @testset "FBDF step!(save_everystep=false) Runtime Allocation Check" for alg in (
            FBDF(), FBDF(time_filter = true), FBDF(time_filter = true, max_order = Val(3)),
        )
        long_prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 100.0))
        integrator = init(
            long_prob, alg, dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6
        )
        # Warm up: take many steps so all caches are initialized and order ramps up
        for _ in 1:50
            step!(integrator)
        end
        # Pre-measurement round to ensure all code paths are fully compiled
        GC.gc()
        for _ in 1:10
            step!(integrator)
        end
        # Actual measurement
        GC.gc()
        allocs = 0
        for _ in 1:10
            allocs += @allocated step!(integrator)
        end
        @test allocs == 0
    end

    @testset "DFBDF step!(save_everystep=false) Runtime Allocation Check" begin
        long_dae_prob = DAEProblem(
            dae_f!, du0, [1.0, 1.0], (0.0, 100.0);
            differential_vars
        )
        integrator = init(
            long_dae_prob, DFBDF(), dt = 0.1, save_everystep = false,
            abstol = 1.0e-6, reltol = 1.0e-6,
            initializealg = BrownFullBasicInit()
        )
        for _ in 1:50
            step!(integrator)
        end
        GC.gc()
        for _ in 1:10
            step!(integrator)
        end
        GC.gc()
        allocs = 0
        for _ in 1:10
            allocs += @allocated step!(integrator)
        end
        @test allocs == 0
    end

    function fbdf_dense_interp_allocs()
        interp_prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 10.0))
        sol = solve(
            interp_prob, FBDF(); abstol = 1.0e-6, reltol = 1.0e-6, dense = true
        )
        out = similar(sol.u[1])
        sol(out, 5.0)
        sol(out, 5.0, Val{1})
        val0 = @allocated sol(out, 5.0)
        val1 = @allocated sol(out, 5.0, Val{1})
        return val0, val1
    end

    @testset "FBDF sol(out, t) Runtime Allocation Check" begin
        @test fbdf_dense_interp_allocs() == (0, 0)
    end

    # Static analysis tests below. These use AllocCheck's check_allocs which
    # analyzes compiled code and may invalidate method caches.
    # Use FullSpecialize to avoid FunctionWrappers dynamic dispatch noise.

    fs_prob = ODEProblem{true, FullSpecialize}(simple_system!, [1.0, 1.0], (0.0, 1.0))
    fs_split_prob = SplitODEProblem(
        SplitFunction(
            ODEFunction{true, FullSpecialize}(f1!),
            ODEFunction{true, FullSpecialize}(f2!)
        ),
        [1.0, 1.0], (0.0, 1.0)
    )
    fs_dae_prob = DAEProblem(
        DAEFunction{true, FullSpecialize}(dae_f!),
        du0, [1.0, 1.0], (0.0, 1.0); differential_vars
    )

    # Test all exported BDF solvers for allocation-free behavior
    # Standard ODE BDF methods
    bdf_solvers = [
        ABDF2(), QNDF1(), QBDF1(), QNDF2(), QBDF2(), QNDF(), QBDF(), FBDF(), MEBDF2(),
    ]

    # IMEX/Split methods need SplitODEProblem
    imex_solvers = [SBDF(order = 2), SBDF2(), SBDF3(), SBDF4(), IMEXEuler(), IMEXEulerARK()]

    # DAE methods need DAEProblem
    dae_solvers = [DABDF2(), DImplicitEuler(), DFBDF()]

    @testset "BDF perform_step! Static Analysis" begin
        for solver in bdf_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    fs_prob, solver, dt = 0.1, save_everystep = false,
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

    @testset "IMEX perform_step! Static Analysis" begin
        for solver in imex_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    fs_split_prob, solver, dt = 0.1, save_everystep = false,
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

    @testset "DAE perform_step! Static Analysis" begin
        for solver in dae_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    fs_dae_prob, solver, dt = 0.1, save_everystep = false,
                    abstol = 1.0e-6, reltol = 1.0e-6,
                    initializealg = BrownFullBasicInit()
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

# A residual the global tracer cannot trace once t >= 0.5 sends a refresh with automatic
# storage through the uncompressed fallback. It must stream columns through reused buffers,
# so its allocation grows linearly with the state count rather than quadratically.
@testset "Uncompressed DAE Jacobian fallback allocates linearly" begin
    function late_branch!(r, du, u, p, t)
        r .= du .+ u
        r[1] += t >= 0.5 && u[1] < 0 ? u[1]^2 : 0.0
        return nothing
    end
    backend = AutoSparse(
        AutoForwardDiff(); sparsity_detector = TracerSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm()
    )
    function refresh_allocation(n)
        prob = DAEProblem(
            late_branch!, -ones(n), ones(n), (0.0, 1.0); differential_vars = trues(n)
        )
        integrator = init(prob, DFBDF(autodiff = backend); initializealg = NoInit())
        integrator.t = 0.75
        nlsolver = integrator.cache.nlsolver
        refresh() = OrdinaryDiffEqDifferentiation.calc_W!(
            nlsolver.cache.W, integrator, nlsolver, integrator.cache, 0.01, false, (true, true)
        )
        refresh()
        refresh()
        GC.gc()
        return @allocated refresh()
    end
    small, large = refresh_allocation(512), refresh_allocation(2048)
    @test large < 6 * small
    @test large < 4_000_000 # a dense 2048 x 2048 matrix alone is 33.5 MB
end
