using OrdinaryDiffEqRKN, Test, RecursiveArrayTools, DiffEqDevTools, Statistics

u0 = fill(0.0, 2)
v0 = ones(2)
function f1_harmonic(dv, v, u, p, t)
    return dv .= -u
end
function f2_harmonic(du, v, u, p, t)
    return du .= v
end
function harmonic_analytic(y0, p, x)
    v0, u0 = y0.x
    return ArrayPartition(-u0 * sin(x) + v0 * cos(x), u0 * cos(x) + v0 * sin(x))
end
ff_harmonic = DynamicalODEFunction(f1_harmonic, f2_harmonic; analytic = harmonic_analytic)
prob = DynamicalODEProblem(ff_harmonic, v0, u0, (0.0, 5.0))

# Methods need BigFloat to test convergence rate
dts = big"1.0" ./ big"2.0" .^ (5:-1:1)
prob_big = DynamicalODEProblem(
    ff_harmonic, [big"1.0", big"1.0"],
    [big"0.0", big"0.0"], (big"0.", big"70.")
)
sim = test_convergence(dts, prob_big, DPRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN6(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 6 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN6FM(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN8(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN12(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 12 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, ERKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, ERKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, ERKN7(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 7 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

# Float64 convergence tests for DPRKN methods
# These tests ensure the CompiledFloats coefficients match the rational coefficients
# (Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/1938)
@testset "Float64 DPRKN convergence" begin
    # Use longer integration time to keep errors above machine precision
    prob_f64 = DynamicalODEProblem(ff_harmonic, ones(2), fill(0.0, 2), (0.0, 20.0))
    dts_f64 = 1.0 ./ 2.0 .^ (2:5)

    sim = test_convergence(dts_f64, prob_f64, DPRKN4())
    @test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1

    sim = test_convergence(dts_f64, prob_f64, DPRKN5())
    @test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1

    sim = test_convergence(dts_f64, prob_f64, DPRKN6())
    @test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1

    sim = test_convergence(dts_f64, prob_f64, DPRKN6FM())
    @test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1

    # DPRKN8 needs larger timesteps and longer integration to avoid hitting machine precision
    prob_f64_long = DynamicalODEProblem(ff_harmonic, ones(2), fill(0.0, 2), (0.0, 50.0))
    dts_f64_large = [1.0, 0.5, 0.25, 0.125]
    sim = test_convergence(dts_f64_large, prob_f64_long, DPRKN8())
    @test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1

    # DPRKN12 is too accurate for Float64 convergence testing (hits machine precision
    # even at dt=1.0), so we only test it with BigFloat (see tests above)
end

sol = solve(prob, Nystrom4(), dt = 1 / 1000)

# Nyström method
dts = 1 .// 2 .^ (9:-1:6)
sim = test_convergence(dts, prob, RKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, Nystrom4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, Nystrom4VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, IRKN3(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 3 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 3 rtol = 1.0e-1
sim = test_convergence(dts, prob, IRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
dts = 1.0 ./ 2.0 .^ (5:-1:0)
sim = test_convergence(dts, prob, Nystrom5VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

# Adaptive methods regression test
sol = solve(prob, FineRKN4())
@test length(sol.u) < 16
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, FineRKN5())
@test length(sol.u) < 14
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN4())
@test length(sol.u) < 25
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN5())
@test length(sol.u) < 38
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN6())
@test length(sol.u) < 20
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN6FM())
@test length(sol.u) < 25
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN8())
@test length(sol.u) < 13
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN12())
@test length(sol.u) < 10
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, ERKN4(), reltol = 1.0e-8)
@test length(sol.u) < 38
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, ERKN5(), reltol = 1.0e-8)
@test length(sol.u) < 34
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, ERKN7(), reltol = 1.0e-8)
@test length(sol.u) < 38
@test SciMLBase.successful_retcode(sol)

u0 = 0.0
v0 = 1.0
function f1_harmonic_nip(v, u, p, t)
    return -u
end
function f2_harmonic_nip(v, u, p, t)
    return v
end

ff_harmonic_nip = DynamicalODEFunction(
    f1_harmonic_nip, f2_harmonic_nip;
    analytic = harmonic_analytic
)
prob = DynamicalODEProblem(ff_harmonic_nip, v0, u0, (0.0, 5.0))

dts = 1 .// 2 .^ (9:-1:6)
sim = test_convergence(dts, prob, RKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, Nystrom4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, Nystrom4VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
@test_broken sim = test_convergence(dts, prob, IRKN3(), dense_errors = true)
@test_broken sim.𝒪est[:l2] ≈ 3 rtol = 1.0e-1
@test_broken sim.𝒪est[:L2] ≈ 3 rtol = 1.0e-1
@test_broken sim = test_convergence(dts, prob, IRKN4(), dense_errors = true)
#@test_broken sim.𝒪est[:l2] ≈ 4 rtol = 1e-1
#@test_broken sim.𝒪est[:L2] ≈ 4 rtol = 1e-1
dts = 1.0 ./ 2.0 .^ (5:-1:0)
sim = test_convergence(dts, prob, Nystrom5VelocityIndependent(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

# Methods need BigFloat to test convergence rate
dts = big"1.0" ./ big"2.0" .^ (5:-1:1)
prob_big = DynamicalODEProblem(
    ff_harmonic_nip, big"1.0", big"0.0",
    (big"0.", big"70.")
)
sim = test_convergence(dts, prob_big, DPRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN6(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 6 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN6FM(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 6 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN8(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 8 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, DPRKN12(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 12 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, ERKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, ERKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob_big, ERKN7(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 7 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

# Adaptive methods regression test
sol = solve(prob, FineRKN4())
@test length(sol.u) < 16
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, FineRKN5())
@test length(sol.u) < 14
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN4())
@test length(sol.u) < 25
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN5())
@test length(sol.u) < 38
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN6())
@test length(sol.u) < 20
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN6FM())
@test length(sol.u) < 25
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN8())
@test length(sol.u) < 13
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, DPRKN12())
@test length(sol.u) < 10
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, ERKN4(), reltol = 1.0e-8)
@test length(sol.u) < 38
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, ERKN5(), reltol = 1.0e-8)
@test length(sol.u) < 34
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, ERKN7(), reltol = 1.0e-8)
@test length(sol.u) < 38
@test SciMLBase.successful_retcode(sol)

# Testing generalized Runge-Kutte-Nyström methods on velocity dependent ODEs with the damped oscillator
println("Out of Place")

# Damped oscillator
prob = ODEProblem(
    DynamicalODEFunction{false}(
        (du, u, p, t) -> -u - 0.5 * du,
        (du, u, p, t) -> du,
        analytic = (du0_u0, p, t) -> ArrayPartition(
            [
                exp(-t / 4) / 15 * (
                    15 * du0_u0[1] * cos(sqrt(15) * t / 4) -
                        sqrt(15) * (du0_u0[1] + 4 * du0_u0[2]) * sin(sqrt(15) * t / 4)
                ),
            ], # du
            [
                exp(-t / 4) / 15 * (
                    15 * du0_u0[2] * cos(sqrt(15) * t / 4) +
                        sqrt(15) * (4 * du0_u0[1] + du0_u0[2]) * sin(sqrt(15) * t / 4)
                ),
            ]
        )
    ),
    ArrayPartition([0.0], [1.0]), # du0, u0
    (0.0, 10.0), # tspan
    SciMLBase.NullParameters(), # p
    SecondOrderODEProblem{false}()
)

dts = 1.0 ./ 2.0 .^ (5:-1:0)
sim = test_convergence(dts, prob, Nystrom4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

# Adaptive methods regression test

sol = solve(prob, FineRKN4())
@test length(sol.u) < 28
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, FineRKN5())
@test length(sol.u) < 20
@test SciMLBase.successful_retcode(sol)

println("In Place")
# Damped oscillator
prob = ODEProblem(
    DynamicalODEFunction{true}(
        (d_du, du, u, p, t) -> @.(d_du = -u - 0.5 * du),
        (d_u, du, u, p, t) -> d_u .= du,
        analytic = (du0_u0, p, t) -> ArrayPartition(
            [
                exp(-t / 4) / 15 * (
                    15 * du0_u0[1] * cos(sqrt(15) * t / 4) -
                        sqrt(15) * (du0_u0[1] + 4 * du0_u0[2]) * sin(sqrt(15) * t / 4)
                ),
            ], # du
            [
                exp(-t / 4) / 15 * (
                    15 * du0_u0[2] * cos(sqrt(15) * t / 4) +
                        sqrt(15) * (4 * du0_u0[1] + du0_u0[2]) * sin(sqrt(15) * t / 4)
                ),
            ]
        )
    ),
    ArrayPartition([0.0], [1.0]), # du0, u0
    (0.0, 10.0), # tspan
    SciMLBase.NullParameters(), # p
    SecondOrderODEProblem{false}()
)

dts = 1.0 ./ 2.0 .^ (5:-1:0)
sim = test_convergence(dts, prob, Nystrom4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN4(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 4 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1
sim = test_convergence(dts, prob, FineRKN5(), dense_errors = true)
@test sim.𝒪est[:l2] ≈ 5 rtol = 1.0e-1
@test sim.𝒪est[:L2] ≈ 4 rtol = 1.0e-1

# Adaptive methods regression test
sol = solve(prob, FineRKN4())
@test length(sol.u) < 28
@test SciMLBase.successful_retcode(sol)
sol = solve(prob, FineRKN5())
@test length(sol.u) < 20
@test SciMLBase.successful_retcode(sol)

# Compare in-place and out-of-place versions
function damped_oscillator(du, u, p, t)
    return -u - 0.5 * du
end
function damped_oscillator!(ddu, du, u, p, t)
    @. ddu = -u - 0.5 * du
    return nothing
end
@testset "in-place vs. out-of-place" begin
    ode_i = SecondOrderODEProblem(
        damped_oscillator!,
        [0.0], [1.0],
        (0.0, 10.0)
    )
    ode_o = SecondOrderODEProblem(
        damped_oscillator,
        [0.0], [1.0],
        (0.0, 10.0)
    )

    @testset "Nystrom4" begin
        alg = Nystrom4()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, dt = dt)
        sol_o = solve(ode_o, alg, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 4 * sol_i.stats.naccept) < 4
    end

    @testset "RKN4" begin
        alg = RKN4()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, dt = dt)
        sol_o = solve(ode_o, alg, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 2 * sol_i.stats.naccept) < 4
    end
    @testset "FineRKN4" begin
        alg = FineRKN4()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 5 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
    end

    @testset "FineRKN5" begin
        alg = FineRKN5()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 7 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        @test_broken sol_i.t ≈ sol_o.t
        @test_broken sol_i.u ≈ sol_o.u
    end

    @testset "DPRKN4" begin
        alg = DPRKN4()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 4 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
    end

    @testset "DPRKN5" begin
        alg = DPRKN5()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 6 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        @test sol_i.t ≈ sol_o.t rtol = 1.0e-5
        @test sol_i.u ≈ sol_o.u rtol = 1.0e-5
    end

    @testset "DPRKN6" begin
        alg = DPRKN6()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test_broken sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 6 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        @test_broken sol_i.t ≈ sol_o.t
        @test_broken sol_i.u ≈ sol_o.u
    end

    @testset "DPRKN6FM" begin
        alg = DPRKN6FM()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 6 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        if VERSION >= v"1.11"
            @test sol_i.t ≈ sol_o.t
            @test sol_i.u ≈ sol_o.u
        else
            @test_broken sol_i.t ≈ sol_o.t
            @test_broken sol_i.u ≈ sol_o.u
        end
    end

    @testset "DPRKN8" begin
        alg = DPRKN8()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 9 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        if VERSION >= v"1.11"
            @test sol_i.t ≈ sol_o.t
            @test sol_i.u ≈ sol_o.u
        else
            @test_broken sol_i.t ≈ sol_o.t
            @test_broken sol_i.u ≈ sol_o.u
        end
    end

    @testset "DPRKN12" begin
        alg = DPRKN12()
        dt = 0.5
        # fixed time step
        sol_i = solve(ode_i, alg, adaptive = false, dt = dt)
        sol_o = solve(ode_o, alg, adaptive = false, dt = dt)
        @test sol_i.t ≈ sol_o.t
        @test sol_i.u ≈ sol_o.u
        @test sol_i.stats.nf == sol_o.stats.nf
        @test sol_i.stats.nf2 == sol_o.stats.nf2
        @test sol_i.stats.naccept == sol_o.stats.naccept
        @test 19 <= sol_i.stats.naccept <= 21
        @test abs(sol_i.stats.nf - 17 * sol_i.stats.naccept) < 4
        # adaptive time step
        sol_i = solve(ode_i, alg)
        sol_o = solve(ode_o, alg)
        @test_broken sol_i.t ≈ sol_o.t
        @test_broken sol_i.u ≈ sol_o.u
    end
end
