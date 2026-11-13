using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    @test sol2.stats.nf <= sol1.stats.nf + 20
end

@testset "Feasible Newton line search" begin
    f!(du, u, p, t) = (du[1] = 1 + u[1] * sqrt(u[1]))
    jac!(J, u, p, t) = (J[1, 1] = 1.5sqrt(u[1]))
    prob = ODEProblem(ODEFunction(f!; jac = jac!), [0.0], (0.0, 0.01))
    sol = solve(
        prob, ImplicitEuler(nlsolve = NLNewton(relax = BackTracking()));
        dt = 0.01, adaptive = false, save_everystep = false
    )
    @test SciMLBase.successful_retcode(sol)
    u = only(sol.u[end])
    @test u > 0
    @test u ≈ 0.01 * (1 + u * sqrt(u)) atol = 1.0e-10 rtol = 0
end
@testset "Newton line search accuracy" begin
    f!(du, u, p, t) = (du[1] = u[1]^2)
    prob = ODEProblem(f!, [1.0], (0.0, 0.0103))
    sol = solve(
        prob, ImplicitEuler(nlsolve = NLNewton(relax = BackTracking()));
        dt = 0.0103, adaptive = false, abstol = 1.0e-12, reltol = 1.0e-12
    )
    reference = setprecision(256) do
        2 / (1 + sqrt(1 - 4big(0.0103)))
    end
    @test SciMLBase.successful_retcode(sol)
    @test only(sol.u[end]) ≈ reference atol = 1.0e-12 rtol = 0
end
