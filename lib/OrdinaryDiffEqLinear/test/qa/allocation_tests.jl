using OrdinaryDiffEqLinear
using OrdinaryDiffEqCore
using SciMLOperators: MatrixOperator
using AllocCheck
using LinearAlgebra
using Test

# Regression bounds for the per-step allocation of the dense (non-Krylov) Magnus
# exponential integrators, which build several NxN stage/commutator matrices every
# step. `linear_perform_step.jl` preallocates those matrices in the cache rather than
# forming new ones on each call; these bounds are well below what master allocates
# (MagnusNC6/MagnusGauss4 at N=50: 1,431,528 / 222,728 bytes/step; at N=200 constant:
# 22,730,984 / 3,524,648 bytes/step) and comfortably above what the preallocated
# version allocates (~60,900 and ~961,900 bytes/step respectively, dominated by the
# scratch matrices ExponentialUtilities.exponential! allocates internally).
@testset "Magnus dense stage-matrix allocation regression" begin
    function magnus_step_alloc(prob, alg; dt = 0.01)
        integrator = init(
            prob, alg; dt = dt, adaptive = false, save_everystep = false
        )
        step!(integrator)
        step!(integrator)
        return @allocated step!(integrator)
    end

    function timedep_update!(A, u, p, t)
        @inbounds for j in axes(A, 2), i in axes(A, 1)
            A[i, j] = (i == j ? -0.1i : 0.001 * sin(t + i * j))
        end
        return nothing
    end

    N = 50
    op = MatrixOperator(zeros(Float64, N, N); update_func! = timedep_update!)
    prob_timedep = ODEProblem(op, ones(Float64, N), (0.0, 0.2))

    Nc = 200
    op_const = MatrixOperator(-0.1 * Matrix{Float64}(I, Nc, Nc))
    prob_const = ODEProblem(op_const, ones(Float64, Nc), (0.0, 0.2))

    @testset "$(typeof(alg)) step allocation bound" for (alg, prob, bound) in (
            (MagnusNC6(), prob_timedep, 150_000),
            (MagnusGauss4(), prob_timedep, 150_000),
            (MagnusNC6(), prob_const, 2_000_000),
            (MagnusGauss4(), prob_const, 2_000_000),
        )
        a = magnus_step_alloc(prob, alg)
        @test a < bound
    end
end

@testset "Linear Allocation Tests" begin
    A = MatrixOperator([-1.0 0.5; 0.0 -2.0])
    prob = ODEProblem(A, [1.0, 1.0], (0.0, 1.0))

    # CayleyEuler excluded: requires matrix-valued state + SplitODEProblem
    linear_solvers = [
        MagnusMidpoint(), LieEuler(), RKMK2(), RKMK4(), LieRK4(),
        CG2(), CG3(), CG4a(), LinearExponential(krylov = :off), MagnusAdapt4(),
    ]

    @testset "Linear perform_step! Static Analysis" begin
        for solver in linear_solvers
            @testset "$(typeof(solver)) perform_step! allocation check" begin
                integrator = init(
                    prob, solver, dt = 0.1, save_everystep = false, adaptive = false
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
