using OrdinaryDiffEqRosenbrock, DiffEqBase, LinearAlgebra, Test

function linear_dae!(du, u, p, t)
    du[1] = -1
    du[2] = u[2] - u[1]
    return nothing
end
linear_dae(u, p, t) = [-1, u[2] - u[1]]

mass_matrix = Diagonal([1.0, 0.0])
callback = ContinuousCallback((u, t, integrator) -> u[1] - 0.3, terminate!)

@testset "Callback-truncated interpolation" begin
    @testset "Mass-matrix DAE" begin
        for f in (linear_dae!, linear_dae), alg in (Rodas4(), Rodas4P(), Rodas5P())
            prob = ODEProblem(ODEFunction(f; mass_matrix), [1.0, 1.0], (0.0, 1.0))
            sol = solve(prob, alg; adaptive = false, dt = 1.0, callback, initializealg = NoInit())

            @test sol.t[end] ≈ 0.7 atol = 1.0e-12
            @test sol.u[end] ≈ [0.3, 0.3] atol = 1.0e-12
            @test sol(0.35) ≈ [0.65, 0.65] atol = 1.0e-12
        end
    end

    @testset "Out-of-place ODE" begin
        prob = ODEProblem((u, p, t) -> u, [1.0], (0.0, 1.0))
        callback = ContinuousCallback((u, t, integrator) -> u[1] - exp(0.7), terminate!)

        for alg in (Rodas4(), Rodas4P(), Rodas5P())
            sol = solve(prob, alg; adaptive = false, dt = 1.0, callback)
            tmiddle = sol.t[end] / 2
            @test sol(tmiddle)[1] ≈ exp(tmiddle) rtol = 1.0e-2
        end
    end
end
