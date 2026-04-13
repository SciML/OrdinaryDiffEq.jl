using OrdinaryDiffEqSSPRK, DiffEqDevTools, Test, Random
import OrdinaryDiffEqLowStorageRK
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

Random.seed!(100)

dts = 1 .// 2 .^ (8:-1:4)
testTol = 0.25

f = (u, p, t) -> cos(t)
prob_ode_sin = ODEProblem(ODEFunction(f; analytic = (u0, p, t) -> sin(t)), 0.0, (0.0, 1.0))

f = (du, u, p, t) -> du[1] = cos(t)
prob_ode_sin_inplace = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> [sin(t)]), [0.0],
    (0.0, 1.0)
)

f = (u, p, t) -> sin(u)
prob_ode_nonlinear = ODEProblem(
    ODEFunction(
        f;
        analytic = (u0, p, t) -> 2 * acot(
            exp(-t) *
                cot(0.5)
        )
    ), 1.0,
    (0.0, 0.5)
)

f = (du, u, p, t) -> du[1] = sin(u[1])
prob_ode_nonlinear_inplace = ODEProblem(
    ODEFunction(
        f;
        analytic = (u0, p, t) -> [
            2 * acot(exp(-t) * cot(0.5)),
        ]
    ),
    [1.0], (0.0, 0.5)
)

test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]

f_ssp = (u, p, t) -> begin
    sin(10t) * u * (1 - u)
end
test_problem_ssp = ODEProblem(f_ssp, 0.1, (0.0, 8.0))
test_problem_ssp_long = ODEProblem(f_ssp, 0.1, (0.0, 1.0e3))

f_ssp_inplace = (du, u, p, t) -> begin
    @. du = sin(10t) * u * (1 - u)
end
test_problem_ssp_inplace = ODEProblem(f_ssp_inplace, rand(3, 3), (0.0, 8.0))

# Test the memory usage, cf. #640
# Note: Basically, the size of the integrator should be the size of the cache
# plus the size of the initial condition, stored is integ.sol.prob.u0.
u0_large = rand(10^6)
prob_ode_large = ODEProblem((du, u, p, t) -> du .= u, u0_large, (0.0, 1.0))

println("SSPRK22")
alg = SSPRK22()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test SSP property of dense output
sol = solve(test_problem_ssp, alg, dt = 1.0)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
sol = solve(test_problem_ssp_inplace, alg, dt = 1.0)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3

println("KYKSSPRK42")
alg = KYKSSPRK42()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)

println("SSPRK33")
alg = SSPRK33()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # This corresponds to Simpson's rule; due to symmetric quadrature nodes,
    # it is of degree 4 instead of 3, as would be expected.
    @test abs(sim.𝒪est[:final] - 1 - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test SSP property of dense output
sol = solve(test_problem_ssp, alg, dt = 1.0)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
sol = solve(test_problem_ssp_inplace, alg, dt = 1.0)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3

println("SSPRK53")
alg = SSPRK53()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4

println("SSPRK53_2N1")
alg = SSPRK53_2N1()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3

# for SSPRK53_2N2 to be in asymptotic range
dts = 1 .// 2 .^ (9:-1:5)
println("SSPRK53_2N2")
alg = SSPRK53_2N2()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 3

dts = 1 .// 2 .^ (9:-1:5)
println("SSPRK53_H")
alg = SSPRK53_H()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = 0.4
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = 0.4
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = 0.4
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 4

#reverting back to original dts
println("SSPRK63")
dts = 1 .// 2 .^ (8:-1:4)
alg = SSPRK63()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)

println("SSPRK73")
alg = SSPRK73()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)

println("SSPRK83")
alg = SSPRK83()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)

println("SSPRK43")
alg = SSPRK43()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test abs(sim.𝒪est[:final] - 1 - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test SSP property of dense output
sol = solve(test_problem_ssp, alg, dt = 8 / 5, adaptive = false)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
sol = solve(test_problem_ssp_inplace, alg, dt = 8 / 5, adaptive = false)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5

println("SSPRK432")
alg = SSPRK432()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # higher order as pure quadrature
    @test abs(sim.𝒪est[:final] - 1 - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test SSP property of dense output
sol = solve(test_problem_ssp, alg, dt = 8 / 5, adaptive = false)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
sol = solve(test_problem_ssp_inplace, alg, dt = 8 / 5, adaptive = false)
@test mapreduce(
    t -> all(0 .<= sol(t) .<= 1), (u, v) -> u && v,
    range(0, stop = 8, length = 50), init = true
)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5

alg = SSPRKMSVS32()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end

println("SSPRKMSVS43")
alg = SSPRKMSVS43()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg) #shows superconvergence to 4th order
    @test abs(sim.𝒪est[:final] - 1 - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end

println("SSPRK932")
alg = SSPRK932()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false, maxiters = 1.0e7
)
@test all(sol.u .>= 0)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5

println("SSPRK54")
alg = SSPRK54()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    # convergence order seems to be worse for this problem
    @test abs(sim.𝒪est[:final] + 0.25 - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    # convergence order seems to be better for this problem
    @test abs(sim.𝒪est[:final] - 0.5 - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)

println("SSPRK104")
alg = SSPRK104()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test sim.𝒪est[:final] ≈ OrdinaryDiffEqSSPRK.alg_order(alg) atol = testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)
# test storage
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 6
integ = init(
    prob_ode_large, alg, dt = 1.0e-2, save_start = false, save_end = false,
    save_everystep = false, alias = ODEAliasSpecifier(alias_u0 = true)
)
@test Base.summarysize(integ) ÷ Base.summarysize(u0_large) <= 5

println("KYK2014DGSSPRK_3S2")
alg = KYK2014DGSSPRK_3S2()
for prob in test_problems_only_time
    sim = test_convergence(dts, prob, alg)
    @test abs(sim.𝒪est[:final] - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
for prob in test_problems_linear
    sim = test_convergence(dts, prob, alg)
    @test abs(sim.𝒪est[:final] - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
for prob in test_problems_nonlinear
    sim = test_convergence(dts, prob, alg)
    @test abs(sim.𝒪est[:final] - OrdinaryDiffEqSSPRK.alg_order(alg)) < testTol
end
# test SSP coefficient
sol = solve(
    test_problem_ssp_long, alg, dt = OrdinaryDiffEqSSPRK.ssp_coefficient(alg),
    dense = false
)
@test all(sol.u .>= 0)

@testset "VectorOfArray/StructArray compatibility" begin
    using RecursiveArrayTools, StaticArrays, StructArrays

    function rhs!(du_voa, u_voa, p, t)
        du = parent(du_voa)
        u = parent(u_voa)
        du .= u
    end

    # StructArray storage
    u = StructArray{SVector{1, Float64}}(ntuple(_ -> [1.0, 2.0], 1))
    ode = ODEProblem(rhs!, VectorOfArray(u), (0, 0.7))
    sol_SA = solve(ode, SSPRK43())

    # Vector{<:SVector} storage
    u = SVector{1, Float64}.([1.0, 2.0])
    ode = ODEProblem(rhs!, VectorOfArray(u), (0, 0.7))
    sol_SV = solve(ode, SSPRK43())

    @test sol_SA ≈ sol_SV
    @test sol_SV.stats.naccept == sol_SA.stats.naccept
end
