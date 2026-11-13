using OrdinaryDiffEqTsit5, ODEProblemLibrary, DiffEqDevTools, Test
using OrdinaryDiffEqRosenbrock, ADTypes
using SciMLBase: successful_retcode

testTol = 0.2
dts = (1 // 2) .^ (7:-1:3)

@testset "Tsit5 Convergence ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (ODEProblemLibrary.prob_ode_linear, ODEProblemLibrary.prob_ode_2Dlinear)[i]
    sim = test_convergence(dts, prob, Tsit5())
    @test abs(sim.𝒪est[:l2] - 5) < testTol + 0.2
end

@testset "AutoTsit5 stiffness switch" begin
    # Robertson (chemical kinetics): stiff, non-oscillatory — endpoint comparison is well-posed
    # (van der Pol limit-cycle phase makes endpoint checks brittle).
    function rober(u, p, t)
        y₁, y₂, y₃ = u
        k₁, k₂, k₃ = p
        return [
            -k₁ * y₁ + k₃ * y₂ * y₃,
            k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃,
            k₂ * y₂^2,
        ]
    end
    prob = ODEProblem(rober, [1.0, 0.0, 0.0], (0.0, 1.0e5), [0.04, 3.0e7, 1.0e4])

    stiffalg = Rosenbrock23(autodiff = AutoFiniteDiff())
    sol = solve(prob, AutoTsit5(stiffalg); maxiters = 10000)
    @test successful_retcode(sol)
    @test length(unique(sol.alg_choice)) > 1  # switched between nonstiff and stiff

    ref = solve(
        prob, Rodas5P(autodiff = AutoFiniteDiff());
        abstol = 1.0e-8, reltol = 1.0e-8,
    )
    @test successful_retcode(ref)
    @test sol.u[end] ≈ ref.u[end] rtol = 1.0e-3 atol = 1.0e-3
end
