using OrdinaryDiffEqExponentialRK
using OrdinaryDiffEqCore
using SciMLBase: FullSpecialize, SplitFunction, ODEFunction
using AllocCheck
using LinearAlgebra
using Test

"""
Allocation tests for OrdinaryDiffEqExponentialRK solvers using AllocCheck.jl.
Tests perform_step! directly (the core stepping function) rather than step!,
since step! includes saving operations that naturally allocate.

Exponential RK solvers operate on split problems: du/dt = L*u + N(u,p,t)
where L is a linear operator and N is the nonlinear part. The problem must
be a SplitODEProblem with MatrixOperator as the linear part.

All exponential RK methods are marked broken=true. These methods compute
matrix exponentials or use Krylov approximations (e.g., Arnoldi iteration),
which inherently allocate. Fixing would require exposing Krylov/Arnoldi
workspace arrays as cache fields.
"""

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
