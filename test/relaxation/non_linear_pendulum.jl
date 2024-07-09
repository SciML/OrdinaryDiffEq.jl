using OrdinaryDiffEq, DiffEqDevTools

printstyled("Non Linear Pendulum\n"; bold = true)

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-sin(u[2]), u[1]]
prob = ODEProblem(
    f,
    [1.0, 0.0],
    (0.0, 1.0))

invariant(x) = x[1]^2/2 -  cos(x[2])

test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)

# Convergence with the method Tsit5()
sim = analyticless_test_convergence(dts, prob, Tsit5(), test_setup)
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ)) 
r = Relaxation(invariant)
sim = analyticless_test_convergence(dts, prob, Tsit5(), test_setup; relaxation = r)
println("order with relaxation with FSAL-R modification: "*string(sim.𝒪est[:final]))