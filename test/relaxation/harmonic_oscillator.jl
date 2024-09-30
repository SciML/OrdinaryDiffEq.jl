using OrdinaryDiffEq, DiffEqDevTools, LinearAlgebra, Test

@testset "Harmonic oscillator" begin

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-u[2],u[1]]
prob = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))

# Convergence with the old method Tsit5()
sim = test_convergence(dts, prob, Tsit5(), adaptive = true)
@test sim.𝒪est[:final] ≈ 5 atol=0.2

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ))
relaxation = Relaxation(norm)
controller = RelaxationController(NonAdaptiveController())
sim_relax = test_convergence(dts, prob, Tsit5(); relaxation, controller)
@test sim_relax.𝒪est[:final] ≈ 5 atol=0.2

end
