using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, SciMLBase, ADTypes, Test

f!(du, u, p, t) = (du .= -u)

@testset "Unresized and empty implicit states" for alg in (Rodas5P, TRBDF2),
        backend in (AutoForwardDiff(), AutoFiniteDiff()), n in (0, 2)
    u0 = ones(n)
    sol = solve(
        ODEProblem{true, SciMLBase.FullSpecialize}(f!, u0, (0.0, 0.25)),
        alg(; autodiff = backend); abstol = 1.0e-10, reltol = 1.0e-10,
        dtmax = 1.0e-4
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol.t[end] == 0.25
    exact = n == 0 ? Float64[] : exp(-0.25) * u0
    @test sol.u[end] ≈ exact rtol = 1.0e-7
end
