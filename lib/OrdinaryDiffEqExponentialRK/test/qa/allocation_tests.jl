using OrdinaryDiffEqExponentialRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, SplitFunction, ODEFunction
using AllocCheck
using LinearAlgebra
using Test

@testset "ExponentialRK Allocation Tests" begin
    A = [-1.0 0.5; 0.0 -2.0]

    # Most exponential RK solvers need SplitODEProblem(MatrixOperator(L), nonlinear!)
    L = MatrixOperator(A)
    function nonlinear!(du, u, p, t)
        du .= 0.0
    end
    f2 = ODEFunction{true, FullSpecialize}(nonlinear!)
    split_fun = SplitFunction(L, f2)
    prob = SplitODEProblem(split_fun, [1.0, 1.0], (0.0, 1.0))

    # Exprb methods require a regular ODEProblem (reject SplitODEProblem)
    function linear_f!(du, u, p, t)
        mul!(du, A, u)
    end
    exprb_prob = ODEProblem{true, FullSpecialize}(linear_f!, [1.0, 1.0], (0.0, 1.0))

    # Solvers that require SplitODEProblem
    exp_solvers = [
        LawsonEuler(), NorsettEuler(), ETD1(), ETDRK2(), ETDRK3(), ETDRK4(),
        HochOst4(), ETD2(),
    ]

    # Exprb methods require regular ODEProblem with full function
    exprb_solvers = [Exprb32(), Exprb43()]

    # Krylov-based methods (Exp4, EPIRK, EXPRB) require additional matrix-free setup
    krylov_solvers = [
        Exp4(), EPIRK4s3A(), EPIRK4s3B(), EPIRK5s3(),
        EXPRB53s3(), EPIRK5P1(), EPIRK5P2(),
    ]

    @testset "ExponentialRK perform_step! Static Analysis" begin
        for solver in exp_solvers
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

        for solver in exprb_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    exprb_prob, solver, dt = 0.1, save_everystep = false,
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

        for solver in krylov_solvers
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

# Runtime allocation guard for the Krylov steppers. Before the column-slice
# updates in perform_step! were `@views`-wrapped, each `X[:, i] .op= ...` read
# materialized a fresh length-`n` column copy, so per-step allocations grew
# linearly with the state size `n`. A view-based (or copy-based) regression can
# only be told apart from the constant-overhead baseline by how it scales with
# `n`, which is version- and platform-robust, unlike an absolute byte ceiling.
# Exp4 uses the symmetric-Jacobian Lanczos path, which keeps a stable Krylov
# subspace, so its per-step allocation is `n`-independent once the slices are
# views; a reintroduced slice copy makes it scale with `n` and trips this test.
@testset "Krylov perform_step! runtime allocations do not scale with state size" begin
    function reaction_diffusion(n)
        dx = 1.0 / (n + 1)
        A = zeros(n, n)
        for i in 1:n
            A[i, i] = -2.0 / dx^2
            i > 1 && (A[i, i - 1] = 1.0 / dx^2)
            i < n && (A[i, i + 1] = 1.0 / dx^2)
        end
        u0 = [sinpi(i * dx) for i in 1:n]
        f! = (du, u, p, t) -> (mul!(du, A, u); @inbounds @. du += u - u^3; nothing)
        function jac!(J, u, p, t)
            copyto!(J, A)
            @inbounds for i in 1:n
                J[i, i] += 1 - 3u[i]^2
            end
            nothing
        end
        ODEProblem(
            ODEFunction{true, FullSpecialize}(f!; jac = jac!, jac_prototype = zeros(n, n)),
            u0, (0.0, 1.0)
        )
    end

    function per_step_allocs(n; nsteps = 40)
        integrator = init(
            reaction_diffusion(n), Exp4(m = 30);
            dt = 1.0e-3, adaptive = false, save_everystep = false, save_start = false
        )
        for _ in 1:5
            step!(integrator)  # warm up compilation and the workspace caches
        end
        return @allocated(
            for _ in 1:nsteps
                step!(integrator)
            end
        ) / nsteps
    end

    small = per_step_allocs(32)
    large = per_step_allocs(256)
    # An unviewed slice copy scales ~8n bytes/step: at n=256 that is several KB
    # over the n=32 baseline, i.e. a large multiple. The view-based code stays
    # essentially flat, so a generous 3x bound cleanly separates the two.
    @test large < 3 * small
end
