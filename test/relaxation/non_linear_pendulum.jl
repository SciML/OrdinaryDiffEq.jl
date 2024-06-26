using OrdinaryDiffEq, DiffEqDevTools

include("relaxation.jl")

printstyled("Non Linear Pendulum\n"; bold = true)

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-sin(u[2]), u[1]]
prob = ODEProblem(
    f,
    [1.0, 0.0],
    (0.0, 1.0))

invariant(x) = x[1]^2/2 -  cos(x[2])

test_setup = Dict(:alg => Vern9(), :reltol => 1e-14, :abstol => 1e-14)

# Convergence with the old method Tsit5()
sim = analyticless_test_convergence(dts, prob, Tsit5(), test_setup)
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))

# Convergence with the new method Tsit5_fors_relaxation() without relaxation
sim = analyticless_test_convergence(dts, prob, Tsit5_for_relaxation(), test_setup)
println("order of convergence of new perform_step! without relaxation: "*string(sim.𝒪est[:final]))

# Convergence with relaxation without FSAL modification, i.e f(uₙ₊₁) ≈ f(uᵧ,ₙ₊₁), before EEst
r = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, invariant))
sim = analyticless_test_convergence(dts, prob, Tsit5_for_relaxation(), test_setup; modif = r)
println("order with relaxation without FSAL modification before EEst: "*string(sim.𝒪est[:final]))

# Convergence with relaxation without FSAL modification, i.e f(uᵧ,ₙ₊₁) ≈ f(uₙ₊₁) , after EEst
r = PerformStepCallback(;postEEst = Relaxation(AlefeldPotraShi, invariant))
sim = analyticless_test_convergence(dts, prob, Tsit5_for_relaxation(), test_setup; modif = r)
println("order with relaxation without FSAL modification after EEst: "*string(sim.𝒪est[:final]))

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ)) 
r = PerformStepCallback(;postEEst = Relaxation(AlefeldPotraShi, invariant, fsal_r))
sim = analyticless_test_convergence(dts, prob, Tsit5_for_relaxation(), test_setup; modif = r)
println("order with relaxation with FSAL-R modification: "*string(sim.𝒪est[:final]))

# Convergence with relaxation with R-FSAL, i.e f(uₙ₊₁) ≈ f(uᵧ,ₙ) + 1/γ ( f(uᵧ,ₙ₊₁) - f(uᵧ,ₙ) )
r = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, invariant, r_fsal), postEEst = fsal_r_int)
sim = analyticless_test_convergence(dts, prob, Tsit5_for_relaxation(), test_setup; modif = r)
println("order with relaxation with R-FSAL modification: "*string(sim.𝒪est[:final]))