using OrdinaryDiffEqTaylorSeries, SparseArrays, Test

function sparse_rhs(u, p, t)
    du = similar(u)
    du[1] = u[3]
    du[2] = 2one(t)
    du[3] = one(t)
    return du
end

function coupled_sparse_rhs(u, p, t)
    du = similar(u)
    du[1] = u[3] + t
    du[2] = 2t - u[1]
    du[3] = 1 + t
    return du
end

@testset "Out-of-place sparse Taylor jets" begin
    u0 = sparsevec([3], [3.0], 3)
    @testset "$(nameof(typeof(solver)))" for solver in (
            ExplicitTaylor(order = Val(5)),
            ExplicitTaylorAdaptiveOrder(min_order = Val(4), max_order = Val(6)),
        )
        @testset "Implicit zeros with nonzero first coefficients" begin
            exact(t) = [3t + t^2 / 2, 2t, 3 + t]
            prob = ODEProblem{false, FullSpecialize}(sparse_rhs, u0, (0.0, 0.1))
            sol = solve(prob, solver; dt = 0.01, adaptive = false)
            @test successful_retcode(sol)
            @test sol.u[end] isa SparseVector
            @test Array(sol.u[end]) ≈ exact(0.1) atol = 1.0e-12 rtol = 0
            @test Array(sol(0.015)) ≈ exact(0.015) atol = 1.0e-12 rtol = 0
        end

        @testset "Coupled polynomial" begin
            exact(t) = [3t + t^2 + t^3 / 6, -t^2 / 2 - t^3 / 3 - t^4 / 24, 3 + t + t^2 / 2]
            prob = ODEProblem{false, FullSpecialize}(coupled_sparse_rhs, u0, (0.0, 0.1))
            sol = solve(prob, solver; dt = 0.005, adaptive = true, abstol = 1.0e-10, reltol = 1.0e-10)
            @test successful_retcode(sol)
            @test axes(sol.u[end]) == axes(u0)
            for t in (0.0, 0.0017, 0.05, 0.1)
                @test Array(sol(t)) ≈ exact(t) atol = 1.0e-11 rtol = 0
            end
        end
    end
end
