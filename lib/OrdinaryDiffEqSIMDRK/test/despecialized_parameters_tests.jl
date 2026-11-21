using OrdinaryDiffEqSIMDRK, Test
using SciMLBase: AutoDespecialize, ODEProblem, solve, successful_retcode

f(u, p, t) = p[1] .* u
prob = ODEProblem{false, AutoDespecialize}(f, [1.0], (0.0, 1.0), [1.0])

for alg in (MER5v2(), MER6v2(), RK6v4())
    sol = solve(prob, alg)
    @test successful_retcode(sol)
    @test only(sol.u[end]) ≈ exp(1) rtol = 1.0e-5
end
