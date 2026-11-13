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

@testset "Local Newton line-search slope" begin
    @testset "forcing $forcing" for forcing in (100.0, 1000.0)
        f!(du, u, p, t) = (du[1] = p - u[1]^3)
        jac!(J, u, p, t) = (J[1, 1] = -3u[1]^2)
        prob = ODEProblem(ODEFunction(f!; jac = jac!), [0.0], (0.0, 1.0), forcing)
        sol = solve(
            prob, ImplicitEuler(
                nlsolve = NLNewton(
                    relax = BackTracking(), max_iter = 100, check_div = false
                )
            );
            dt = 1.0, adaptive = false, abstol = 1.0e-10, reltol = 1.0e-10
        )
        reference = setprecision(256) do
            u = cbrt(big(forcing))
            for _ in 1:20
                u -= (u + u^3 - forcing) / (1 + 3u^2)
            end
            u
        end
        @test SciMLBase.successful_retcode(sol)
        @test only(sol.u[end]) ≈ reference atol = 1.0e-8 rtol = 0
    end
end

@testset "Large finite Newton line search" begin
    f!(du, u, p, t) = (du[1] = -u[1])
    jac!(J, u, p, t) = (J[1, 1] = -1)
    prob = ODEProblem(ODEFunction(f!; jac = jac!), [1.0e200], (0.0, 0.1))
    sol = solve(
        prob, ImplicitEuler(nlsolve = NLNewton(relax = BackTracking()));
        dt = 0.1, adaptive = false, abstol = 1.0e194, reltol = 1.0e-6
    )
    reference = big(1.0e200) / (1 + big(0.1))
    @test SciMLBase.successful_retcode(sol)
    @test only(sol.u[end]) ≈ reference rtol = 1.0e-12
end
