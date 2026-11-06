using OrdinaryDiffEqDefault, Test, LinearAlgebra, SparseArrays

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
            prototypes = (
                prototypes..., (Symmetric(ones(n, n)), 6),
                (Hermitian(ones(n, n)), 6),
                (Diagonal(view(ones(n), :)), 5),
            )
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
