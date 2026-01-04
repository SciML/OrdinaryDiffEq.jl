using OrdinaryDiffEqLowOrderRK
f_ssp = (u, p, t) -> begin
    sin(10t) * u * (1 - u)
end
test_problem_ssp = ODEProblem(f_ssp, 0.1, (0.0, 8.0))
test_problem_ssp_long = ODEProblem(f_ssp, 0.1, (0.0, 1.0e3))

# test SSP coefficient for explicit Euler
alg = Euler()
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqLowOrderRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqLowOrderRK.ssp_coefficient(alg) + 1.0e-3,
    dense = false
)
@test any(sol.u .< 0)
