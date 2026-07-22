using OrdinaryDiffEqNonlinearSolve: HomotopyNonlinearSolveAlg
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF
using SciMLBase
using NonlinearSolve
using ADTypes: AutoFiniteDiff
using StaticArrays
using LinearAlgebra
using Test

function vdp!(du, u, p, t)
    du[1] = u[2]
    du[2] = p[1] * ((1 - u[1]^2) * u[2] - u[1])
    return nothing
end
vdp(u, p, t) = SA[u[2], p[1] * ((1 - u[1]^2) * u[2] - u[1])]

function cached_homotopy_step_allocations!(integrator)
    step!(integrator)
    step!(integrator)
    GC.gc()
    return @allocated step!(integrator)
end

@testset "stiff van der Pol, DIRK and multistep forms" begin
    prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, 0.3), [1.0e3])
    ref = solve(prob, TRBDF2(), abstol = 1.0e-12, reltol = 1.0e-12)
    hnls = HomotopyNonlinearSolveAlg()

    for alg in (
            ImplicitEuler(nlsolve = hnls), TRBDF2(nlsolve = hnls),
            Trapezoid(nlsolve = hnls), QNDF(nlsolve = hnls), FBDF(nlsolve = hnls),
        )
        sol = solve(prob, alg, abstol = 1.0e-6, reltol = 1.0e-6)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test maximum(abs.(sol.u[end] .- ref.u[end])) < 1.0e-2
    end
end

@testset "continuation cache is reused across stages" begin
    prob = ODEProblem((du, u, p, t) -> (du[1] = t - u[1]; nothing), [1.0], (0.0, 1.0))
    nlsolve_algs = (
        HomotopyNonlinearSolveAlg(),
        HomotopyNonlinearSolveAlg(
            KantorovichHomotopy(
                inner = NewtonRaphson(autodiff = AutoFiniteDiff()), strict = false
            )
        ),
    )
    for nlsolve_alg in nlsolve_algs
        integrator = init(
            prob, ImplicitEuler(nlsolve = nlsolve_alg);
            adaptive = false, dt = 0.1, save_everystep = false
        )
        continuation_cache = integrator.cache.nlsolver.cache.continuation_cache
        @test continuation_cache !== nothing

        step!(integrator)
        @test integrator.cache.nlsolver.cache.continuation_cache === continuation_cache
        step!(integrator)
        @test integrator.cache.nlsolver.cache.continuation_cache === continuation_cache
        first_step = (1.0 + 0.1 * 0.1) / 1.1
        @test integrator.u[1] ≈ (first_step + 0.1 * 0.2) / 1.1
        VERSION >= v"1.11" && @test cached_homotopy_step_allocations!(integrator) == 0
    end
end

@testset "continuation cache is rebuilt after integrator resize" begin
    fresize! = function (du, u, p, t)
        for i in eachindex(u)
            du[i] = -u[i]
        end
        return nothing
    end
    prob = ODEProblem(fresize!, [1.0], (0.0, 0.3))
    integrator = init(
        prob, ImplicitEuler(nlsolve = HomotopyNonlinearSolveAlg());
        adaptive = false, dt = 0.1
    )
    step!(integrator)
    old_cache = integrator.cache.nlsolver.cache.continuation_cache

    resize!(integrator, 2)
    integrator.u .= 1.0
    step!(integrator)

    @test length(integrator.u) == 2
    @test integrator.cache.nlsolver.cache.continuation_cache !== old_cache
    @test !integrator.cache.nlsolver.cache.needs_rebuild
end

@testset "out-of-place with StaticArrays" begin
    prob = ODEProblem(vdp, SA[2.0, 0.0], (0.0, 0.3), SA[1.0e3])
    ref = solve(prob, TRBDF2(), abstol = 1.0e-12, reltol = 1.0e-12)
    hnls = HomotopyNonlinearSolveAlg()
    for alg in (ImplicitEuler(nlsolve = hnls), TRBDF2(nlsolve = hnls))
        sol = solve(prob, alg, abstol = 1.0e-6, reltol = 1.0e-6)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test maximum(abs.(sol.u[end] .- ref.u[end])) < 1.0e-2
    end
end

@testset "ArcLengthContinuation as the continuation solver" begin
    prob = ODEProblem(vdp!, [2.0, 0.0], (0.0, 0.3), [1.0e3])
    ref = solve(prob, TRBDF2(), abstol = 1.0e-12, reltol = 1.0e-12)
    arc = HomotopyNonlinearSolveAlg(
        ArcLengthContinuation(inner = NewtonRaphson(autodiff = AutoFiniteDiff()))
    )
    sol = solve(prob, ImplicitEuler(nlsolve = arc), abstol = 1.0e-6, reltol = 1.0e-6)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test maximum(abs.(sol.u[end] .- ref.u[end])) < 1.0e-2
end

@testset "Robertson" begin
    function rober!(du, u, p, t)
        y1, y2, y3 = u
        k1, k2, k3 = p
        du[1] = -k1 * y1 + k3 * y2 * y3
        du[2] = k1 * y1 - k2 * y2^2 - k3 * y2 * y3
        du[3] = k2 * y2^2
        return nothing
    end
    prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e5), [0.04, 3.0e7, 1.0e4])
    ref = solve(prob, FBDF(), abstol = 1.0e-12, reltol = 1.0e-12)
    hnls = HomotopyNonlinearSolveAlg()
    for alg in (TRBDF2(nlsolve = hnls), FBDF(nlsolve = hnls))
        sol = solve(prob, alg, abstol = 1.0e-8, reltol = 1.0e-8)
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test maximum(abs.(sol.u[end] .- ref.u[end])) < 1.0e-5
    end
end

@testset "nonsingular mass matrix" begin
    mmf = ODEFunction(
        (du, u, p, t) -> (du[1] = -u[1] + u[2]; du[2] = -10u[2]; nothing);
        mass_matrix = Diagonal([2.0, 3.0])
    )
    prob = ODEProblem(mmf, [1.0, 1.0], (0.0, 1.0))
    sol = solve(
        prob, ImplicitEuler(nlsolve = HomotopyNonlinearSolveAlg()),
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    refsol = solve(prob, ImplicitEuler(), abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test maximum(abs.(sol.u[end] .- refsol.u[end])) < 1.0e-8
end

@testset "blow-up problem passes through fold via step rejection" begin
    # u' = u², u(0) = 1 blows up at t = 1; the implicit Euler stage equation has a
    # fold at effective step size 1/(4u), so large-dt attempts must be rejected by
    # the continuation (no consistent solution exists) and retried smaller
    prob = ODEProblem((u, p, t) -> u .^ 2, [1.0], (0.0, 0.9))
    sol = solve(
        prob, ImplicitEuler(nlsolve = HomotopyNonlinearSolveAlg()),
        abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test abs(sol.u[end][1] - 1 / (1 - 0.9)) < 0.05
end

@testset "DAE problems are rejected" begin
    dae = DAEProblem((du, u, p, t) -> du .- u, [1.0], [1.0], (0.0, 1.0))
    @test_throws ArgumentError solve(
        dae, DABDF2(nlsolve = HomotopyNonlinearSolveAlg())
    )
end
