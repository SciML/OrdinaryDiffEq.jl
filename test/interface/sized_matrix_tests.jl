using OrdinaryDiffEq, StaticArrays, Test

function dudt!(du, σ, p, t)
    du .= -σ
end

η0 = 1
τ = 1
α = 0.8

## Exercise code paths that are in-place but not `<:Array` but not AbstractVector
σ0 = SizedMatrix{3, 3}([1.0 0.0 0.0; 0.0 2.0 0.0; 3.0 0.0 0.0])

prob_giesekus = ODEProblem(dudt!, σ0, (0.0, 2.0))
prob_normal = ODEProblem(dudt!, Matrix(σ0), (0.0, 2.0))


if VERSION >= v"1.9"
    @testset "giesekus" begin
    @testset "$alg" for alg in [
        Rosenbrock23(),
        Rodas4(),
        Rodas4P(),
        Rodas5(),
        Rodas5P(),
        Tsit5(),
        Vern6(),
        Vern7(),
        Vern8(),
        Vern9(),
        DP5()
    ]
        sol = solve(prob_giesekus, alg, saveat = 0.2, abstol = 1e-10, reltol = 1e-10)
        sol_normal = solve(prob_giesekus, alg, saveat = 0.2, abstol = 1e-10, reltol = 1e-10)
        @test sol.retcode == sol_normal.retcode == ReturnCode.Success
        @test Array(sol) ≈ Array(sol_normal)
        # give ourselves some stats wiggle room
        @test sol.stats.nreject <= 2*sol_normal.stats.nreject + 1
        @test sol.stats.naccept <= 2*sol_normal.stats.naccept + 1
    end
    end
end
