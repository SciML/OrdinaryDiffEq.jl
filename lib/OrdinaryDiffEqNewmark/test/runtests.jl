using OrdinaryDiffEqNewmark, Test, RecursiveArrayTools, DiffEqDevTools, Statistics

# Newmark methods with harmonic oscillator
@testset "Harmonic Oscillator" begin
    u0 = fill(0.0, 2)
    v0 = ones(2)
    function f1_harmonic!(dv, v, u, p, t)
        dv .= -u
    end
    function f2_harmonic!(du, v, u, p, t)
        du .= v
    end
    function harmonic_analytic(y0, p, x)
        v0, u0 = y0.x
        ArrayPartition(-u0 * sin(x) + v0 * cos(x), u0 * cos(x) + v0 * sin(x))
    end

    ff_harmonic! = DynamicalODEFunction(f1_harmonic!, f2_harmonic!; analytic = harmonic_analytic)
    prob = DynamicalODEProblem(ff_harmonic!, v0, u0, (0.0, 5.0))
    dts = 1.0 ./ 2.0 .^ (5:-1:0)

    sim = test_convergence(dts, prob, Newmark(), dense_errors = true)
    @test sim.𝒪est[:l2]≈2 rtol=1e-1
end

# Newmark methods with damped oscillator
@testset "Damped Oscillator" begin
    # function damped_oscillator(du, u, p, t)
    #     return -u - 0.5 * du
    # end
    # function damped_oscillator!(ddu, du, u, p, t)
    #     @. ddu = -u - 0.5 * du
    #     return nothing
    # end
    function damped_oscillator_analytic(du0_u0, p, t)
        OrdinaryDiffEq.SciMLBase.ArrayPartition(
            [
                exp(-t / 4) / 15 * (15 * du0_u0[1] * cos(sqrt(15) * t / 4) -
                sqrt(15) * (du0_u0[1] + 4 * du0_u0[2]) * sin(sqrt(15) * t / 4))
            ], # du
            [
                exp(-t / 4) / 15 * (15 * du0_u0[2] * cos(sqrt(15) * t / 4) +
                sqrt(15) * (4 * du0_u0[1] + du0_u0[2]) * sin(sqrt(15) * t / 4))
            ]
    end
    ff_harmonic_damped! = DynamicalODEFunction((ddu, v, u, p, t) -> ddu = -u - 0.5 * v,
        (du, v, u, p, t) -> du = v,
        analytic = (du0_u0, p, t) -> damped_oscillator_analytic
        )
    )

    prob = DynamicalODEProblem(ff_harmonic_damped!, [0.0], [1.0], (0.0, 10.0))
    dts = 1.0 ./ 2.0 .^ (5:-1:0)

    sim = test_convergence(dts, prob, Newmark(), dense_errors = true)
    @test sim.𝒪est[:l2]≈2 rtol=1e-1
end

# @testset "in-place vs. out-of-place" begin
#     ode_i = SecondOrderODEProblem(damped_oscillator!,
#         [0.0], [1.0],
#         (0.0, 10.0))
#     ode_o = SecondOrderODEProblem(damped_oscillator,
#         [0.0], [1.0],
#         (0.0, 10.0))

#     @testset "NewmarkBeta" begin
#         alg = NewmarkBeta()
#         dt = 0.5
#         # fixed time step
#         sol_i = solve(ode_i, alg, dt = dt)
#         sol_o = solve(ode_o, alg, dt = dt)
#         @test sol_i.t ≈ sol_o.t
#         @test sol_i.u ≈ sol_o.u
#         @test sol_i.destats.nf == sol_o.destats.nf
#         @test sol_i.destats.nf2 == sol_o.destats.nf2
#         @test sol_i.destats.naccept == sol_o.destats.naccept
#         @test 19 <= sol_i.destats.naccept <= 21
#         @test abs(sol_i.destats.nf - 4 * sol_i.destats.naccept) < 4
#     end
# end
