using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

@testset "Fixed-step Newton does not stop at the first iterate" begin
    f!(du, u, p, t) = (du[1] = -u[1]^3 + 10 * sin(5t); du[2] = -100 * (u[2] - u[1]^2); nothing)
    sol = solve(ODEProblem(f!, [1.0, 1.0], (0.0, 5.0)), ImplicitEuler(); dt = 0.02, adaptive = false)
    du = similar(sol.u[1])
    residual = maximum(2:length(sol.t)) do n
        f!(du, sol.u[n], nothing, sol.t[n])
        sqrt(sum(abs2, sol.u[n] - sol.u[n - 1] - (sol.t[n] - sol.t[n - 1]) * du) / sum(abs2, sol.u[n]))
    end
    @test residual < 1.0e-4
end

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
