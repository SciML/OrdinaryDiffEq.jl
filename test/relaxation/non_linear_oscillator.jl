using OrdinaryDiffEq, DiffEqDevTools, LinearAlgebra, Test

@testset "Non linear harmonic oscillator" begin

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-u[2]/(u[1]^2 + u[2]^2),u[1]/(u[1]^2 + u[2]^2)]
prob = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))


# Convergence with the method Tsit5()
sim = test_convergence(dts, prob, Tsit5())
@test sim.𝒪est[:final] ≈ 5.4 atol=0.2

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ))
relaxation = Relaxation(norm)
controller = RelaxationController(NonAdaptiveController())
sim_relax = test_convergence(dts, prob, Tsit5(); relaxation, controller, adaptive=true)
@test sim.𝒪est[:final] ≈ 5.4 atol=0.2

end
