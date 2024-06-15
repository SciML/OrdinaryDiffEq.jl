using OrdinaryDiffEq, DiffEqDevTools

include("relaxation.jl")

printstyled("Harmonic Oscillator\n"; bold = true)

dts = (1 / 2) .^ (6:-1:4)

f = (u, p, t) -> [-u[2],u[1]]
prob = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> [cos(t), sin(t)]),
    [1.0, 0.0],
    (0.0, 1.0))

# Convergence with the old method Tsit5()
sim = test_convergence(dts, prob, Tsit5())
println("order of convergence of older perform_step! : "*string(sim.𝒪est[:final]))

# Convergence with the new method Tsit5_fors_relaxation() without relaxation
sim = test_convergence(dts, prob, Tsit5_for_relaxation())
println("order of convergence of new perform_step! without relaxation: "*string(sim.𝒪est[:final]))

# Convergence with relaxation without FSAL modification, i.e f(uₙ₊₁) ≈ f(uᵧ,ₙ₊₁), before EEst
r = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> norm(x)))
sim = test_convergence(dts, prob, Tsit5_for_relaxation(); modif = r)
println("order with relaxation without FSAL modification before EEst: "*string(sim.𝒪est[:final]))

# Convergence with relaxation without FSAL modification, i.e f(uᵧ,ₙ₊₁) ≈ f(uₙ₊₁) , after EEst
r = PerformStepCallback(;postEEst = Relaxation(AlefeldPotraShi, x-> norm(x)))
sim = test_convergence(dts, prob, Tsit5_for_relaxation(); modif = r)
println("order with relaxation without FSAL modification after EEst: "*string(sim.𝒪est[:final]))

# Convergence with relaxation with FSAL-R, i.e  f(uᵧ,ₙ₊₁) ≈ f(uᵧ,ₙ) + γ ( f(uₙ₊₁) - f(uᵧ,ₙ)) 
r = PerformStepCallback(;postEEst = Relaxation(AlefeldPotraShi, x-> norm(x), fsal_r))
sim = test_convergence(dts, prob, Tsit5_for_relaxation(); modif = r)
println("order with relaxation with FSAL-R modification: "*string(sim.𝒪est[:final]))

# Convergence with relaxation with R-FSAL, i.e f(uₙ₊₁) ≈ f(uᵧ,ₙ) + 1/γ ( f(uᵧ,ₙ₊₁) - f(uᵧ,ₙ) )
r = PerformStepCallback(;poststep = Relaxation(AlefeldPotraShi, x-> norm(x), r_fsal), postEEst = fsal_r_int)
sim = test_convergence(dts, prob, Tsit5_for_relaxation(); modif = r)
println("order with relaxation with R-FSAL modification: "*string(sim.𝒪est[:final]))