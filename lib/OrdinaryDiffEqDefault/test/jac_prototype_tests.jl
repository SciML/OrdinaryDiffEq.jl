using OrdinaryDiffEqDefault, Test, LinearAlgebra, SparseArrays
using BandedMatrices: BandedMatrix
using BlockBandedMatrices: BandedBlockBandedMatrix, BlockBandedMatrix
using BlockDiagonals: BlockDiagonal
using FastAlmostBandedMatrices: AlmostBandedMatrix
using LinearSolve: DiagonalFactorization, GenericFactorization, KrylovJL_GMRES, LUFactorization

@testset "Default stiff Jacobian storage" begin
    for n in (500, 501)
        rates = collect(range(1.0, 1000.0; length = n))
        decay!(du, u, p, t) = (@. du = -p * u)
        prototypes = (
            (Diagonal(ones(n)), 5),
            (Tridiagonal(ones(n - 1), ones(n), ones(n - 1)), 5),
            (sparse(Diagonal(ones(n))), 5),
            (ones(n, n), n > 500 ? 6 : 5),
            (nothing, n > 500 ? 6 : 5),
        )
        if n > 500
            blocks = [fill(2, (n - 1) ÷ 2); 1]
            block_pattern = Tridiagonal(ones(n - 1), ones(n), ones(n - 1))
            prototypes = (
                prototypes..., (Symmetric(ones(n, n)), 6),
                (Hermitian(ones(n, n)), 6),
                (SymTridiagonal(ones(n), ones(n - 1)), 6),
                (Diagonal(view(ones(n), :)), 5),
                (Bidiagonal(ones(n), ones(n - 1), :U), 5),
                (BandedMatrix(ones(n, n), (2, 2)), 5),
                (BlockBandedMatrix(block_pattern, blocks, blocks, (1, 1)), 5),
                (BandedBlockBandedMatrix(block_pattern, blocks, blocks, (1, 1), (1, 1)), 5),
                (BlockDiagonal([ones(b, b) for b in blocks]), 5),
                (AlmostBandedMatrix(BandedMatrix(ones(n, n), (2, 2)), ones(1, n)), 5),
            )
            if isdefined(SparseArrays, :fixed)
                prototypes = (prototypes..., (SparseArrays.fixed(sparse(Diagonal(ones(n)))), 6))
            end
        end
        for (jac_prototype, slot) in prototypes
            @testset "n=$n, $(typeof(jac_prototype))" begin
                prob = ODEProblem(
                    ODEFunction(decay!; jac_prototype), ones(n), (0.0, 0.01), rates
                )
                sol = solve(
                    prob, DefaultImplicitODEAlgorithm(); abstol = 1.0e-9, reltol = 1.0e-7
                )
                @test successful_retcode(sol)
                @test all(==(slot), sol.alg_choice)
                @test maximum(abs, sol.u[end] .- exp.(-rates .* 0.01)) < 1.0e-5
            end
        end
    end

    n = 501
    rates = collect(range(1.0, 1000.0; length = n))
    decay_mass!(du, u, p, t) = (@. du = -p * u)
    prob = ODEProblem(
        ODEFunction(decay_mass!; jac_prototype = Diagonal(ones(n)), mass_matrix = 2I),
        ones(n), (0.0, 0.01), rates
    )
    sol = solve(prob; abstol = 1.0e-9, reltol = 1.0e-7)
    @test successful_retcode(sol)
    @test all(==(5), sol.alg_choice)
    @test maximum(abs, sol.u[end] .- exp.(-rates .* 0.01 ./ 2)) < 1.0e-5
end

@testset "Explicit linsolve on both large-system paths" begin
    n = 501
    rates = collect(range(1.0, 1000.0; length = n))
    decay!(du, u, p, t) = (@. du = -p * u)
    band_qr(A) = qr(A)
    for (prototype, slot, direct) in (
            (Diagonal(ones(n)), 5, DiagonalFactorization()),
            (BandedMatrix(ones(n, n), (2, 2)), 5, GenericFactorization(; fact_alg = band_qr)),
            (ones(n, n), 6, LUFactorization()),
        )
        f = ODEFunction(decay!; jac_prototype = prototype, mass_matrix = 2I)
        prob = ODEProblem(f, ones(n), (0.0, 0.01), rates)
        for linsolve in (direct, KrylovJL_GMRES()), alg in (
                    DefaultODEAlgorithm(; linsolve, stiffalgfirst = true),
                    DefaultImplicitODEAlgorithm(; linsolve),
                )
            integrator = init(prob, alg; abstol = 1.0e-9, reltol = 1.0e-7)
            @test all(member -> member.linsolve === linsolve, integrator.alg.algs[3:end])
            sol = solve!(integrator)
            @test successful_retcode(sol)
            @test all(==(slot), sol.alg_choice)
            @test (sol.stats.njacs > 0) == (linsolve === direct)
            @test maximum(abs, sol.u[end] .- exp.(-rates .* 0.01 ./ 2)) < 1.0e-5
        end
    end
end
