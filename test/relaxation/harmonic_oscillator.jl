using OrdinaryDiffEq, DiffEqDevTools

printstyled("Harmonic Oscillator\n"; bold = true)

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-u[2],u[1]]
prob = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))

invariant(x) = norm(x)

# Convergence with the old method Tsit5()
sim = test_convergence(dts, prob, Tsit5())
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ)) 
r = Relaxation(invariant)
sim = test_convergence(dts, prob, Tsit5(); relaxation = r)
println("order with relaxation with FSAL-R modification: "*string(sim.𝒪est[:final]))

