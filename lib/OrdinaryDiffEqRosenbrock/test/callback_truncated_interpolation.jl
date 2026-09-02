using OrdinaryDiffEqRosenbrock, DiffEqBase, LinearAlgebra, Test

function linear_dae!(du, u, p, t)
    du[1] = -u[1]
    du[2] = u[2] + u[1] - 1
    return nothing
end
linear_dae(u, p, t) = [-u[1], u[2] + u[1] -1]

mass_matrix = Diagonal([1.0, 0.0])
callback = ContinuousCallback((u, t, integrator) -> t - 0.7, terminate!)

@testset "Callback-truncated interpolation" begin
    @testset "Mass-matrix DAE" begin
        for f in (linear_dae!, linear_dae), alg in (Rodas4(), Rodas4P(), Rodas5P())
            prob = ODEProblem(ODEFunction(f; mass_matrix), [1.0, 0.0], (0.0, 1.0))
            sol = solve(prob, alg; callback, initializealg = CheckInit())
            @test sol.retcode == SciMLBase.ReturnCode.Terminated

            @test sol.t[end] ≈ 0.7 atol = 1.0e-12
            @test sol.u[end] ≈ [exp(-0.7), 1-exp(-0.7)] atol = 1.0e-4
            @testset "Dense interpolation" begin
                t_dense = range(0, stop=sol.t[end], length=101)
                for ti in t_dense
                    @test sol(ti) ≈ [exp.(-ti), 1 .- exp.(-ti)] atol = 1.0e-2
                end
            end
        end
    end

    @testset "Out-of-place ODE" begin
        prob = ODEProblem((u, p, t) -> u, [1.0], (0.0, 1.0))
        for alg in (Rodas4(), Rodas4P(), Rodas5P())
            sol = solve(prob, alg; callback)
            @test sol.retcode == SciMLBase.ReturnCode.Terminated
            tmiddle = sol.t[end] / 2
            @test sol(tmiddle)[1] ≈ exp(tmiddle) rtol = 1.0e-2
            @testset "Dense interpolation" begin
                t_dense = range(0, stop=sol.t[end], length=101)
                for ti in t_dense
                    @test sol(ti) ≈ [exp.(ti)] atol = 1.0e-2
                end
            end
        end
    end
end