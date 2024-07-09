using OrdinaryDiffEq, DiffEqDevTools

printstyled("Time-dependent harmonic oscillator\n"; bold = true)

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-(1+0.5 * sin(t))*u[2], (1+0.5 * sin(t))*u[1]]
prob = ODEProblem(
    ODEFunction(f; 
                analytic = (u0, p, t)->[cos(0.5)*cos(t-0.5*cos(t))-sin(0.5)*sin(t-0.5*cos(t)), 
                                        sin(0.5)*cos(t-0.5*cos(t))+cos(0.5)*sin(t-0.5*cos(t))]),
    [1.0, 0.0],
    (0.0, 1.0))

invariant(x) = norm(x)

# Convergence with the method Tsit5()
sim = test_convergence(dts, prob, Tsit5(), adaptative = true)
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ)) 
r = Relaxation(invariant)
sim = test_convergence(dts, prob, Tsit5(); relaxation = r, adaptative = true)
println("order with relaxation with FSAL-R modification: "*string(sim.𝒪est[:final]))