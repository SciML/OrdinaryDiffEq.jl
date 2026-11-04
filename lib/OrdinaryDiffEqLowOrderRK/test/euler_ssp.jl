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

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/2189
@testset "Fractional tstops with integer tspan" begin
    tstops = [0.33, 0.66, 1.0]
    sol = solve(ODEProblem((du, u, p, t) -> du .= u, [1.0], (0, 1)), Euler(); tstops)
    @test OrdinaryDiffEqLowOrderRK.SciMLBase.successful_retcode(sol)
    @test sol.t == [0.0, 0.33, 0.66, 1.0]
    @test sol.u[end] ≈ [2.370326]
    @test eltype(sol.t) === Float64
    for stops in ((0.33, 0.66, 1.0), Any[0, 0.33, 0.66, 1])
        sol = solve(ODEProblem((du, u, p, t) -> du .= u, [1.0], (0, 1)), Euler(); tstops = stops)
        @test sol.t == [0.0, 0.33, 0.66, 1.0]
        @test sol.u[end] ≈ [2.370326]
    end
    scalar_stop = solve(ODEProblem((du, u, p, t) -> du .= u, [1.0], (0, 1)), Euler(); tstops = 0.33)
    @test OrdinaryDiffEqLowOrderRK.SciMLBase.successful_retcode(scalar_stop)
    @test scalar_stop.t == [0.0, 0.33, 1.0]
    @test scalar_stop.u[end] ≈ [2.2211]
    problem_stops = ODEProblem((du, u, p, t) -> du .= u, [1.0], (0, 1); tstops)
    problem_sol = solve(problem_stops, Euler())
    @test problem_sol.t == [0.0, 0.33, 0.66, 1.0]
    @test problem_sol.u[end] ≈ [2.370326]
    @test OrdinaryDiffEqLowOrderRK.SciMLBase.successful_retcode(problem_sol)
    @test eltype(problem_sol.t) === Float64
    overridden_sol = solve(problem_stops, Euler(); tstops = [1])
    @test overridden_sol.t == [0, 1]
    @test overridden_sol.u[end] == [2.0]
    @test eltype(overridden_sol.t) === Int
    for stops in ([1 // 3, 2 // 3, 1 // 1], [1, 2])
        final_time = last(stops)
        sol = solve(ODEProblem((du, u, p, t) -> du .= u, [1.0], (0, final_time)), Euler(); tstops = stops)
        @test OrdinaryDiffEqLowOrderRK.SciMLBase.successful_retcode(sol)
        @test eltype(sol.t) === eltype(stops)
        @test sol.t == [0; stops]
        @test sol.u[end] ≈ [prod(1 .+ diff([0; stops]))]
    end
    large_start = 2^53
    for stops in ([0.5], [Float64(large_start)], [Float64(large_start + 2)])
        exact = solve(
            ODEProblem((du, u, p, t) -> du .= u, [1.0], (large_start, large_start + 1)),
            Euler(); dt = 1, tstops = stops
        )
        @test exact.t == [large_start, large_start + 1]
        @test exact.u[end] == [2.0]
        @test eltype(exact.t) === Int
    end
    integral_stop = solve(
        ODEProblem((du, u, p, t) -> du .= u, [1.0], (large_start, large_start + 3)),
        Euler(); dt = 2, tstops = [Float64(large_start + 2)]
    )
    @test integral_stop.t == [large_start, large_start + 2, large_start + 3]
    @test integral_stop.u[end] == [6.0]
    @test eltype(integral_stop.t) === Int
    backward = solve(ODEProblem((du, u, p, t) -> du .= u, [1.0], (1, 0)), Euler(); tstops = [0.66, 0.33, 0.0])
    @test backward.t == [1.0, 0.66, 0.33, 0.0]
    @test backward.u[end] ≈ [0.66 * 0.67 * 0.67]
end
