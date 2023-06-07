using OrdinaryDiffEq, DiffEqDevTools, DiffEqBase, Test

import ODEProblemLibrary: prob_ode_bigfloatlinear,
    prob_ode_linear,
    prob_ode_2Dlinear,
    prob_ode_bigfloat2Dlinear

probbig = prob_ode_bigfloat2Dlinear
probnum = prob_ode_linear
probnumbig = prob_ode_bigfloatlinear
prob = prob_ode_2Dlinear

dts = (1 / 2) .^ (7:-1:4)
testTol = 0.2
bools = Vector{Bool}(undef, 0)

### BS3()
println("BS3")
sim = test_convergence(dts, probnum, BS3())
@test abs.(sim.𝒪est[:l2] - 3) < testTol
sim = test_convergence(dts, prob, BS3())
@test abs.(sim.𝒪est[:l2] - 3) < testTol

tabalg = ExplicitRK(tableau = constructBogakiShampine3())
sol1 = solve(probnum, BS3(), dt = 1 / 2^1, adaptive = false, save_everystep = false)
sol2 = solve(probnum, tabalg, dt = 1 / 2^1, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(prob, BS3(), dt = 1 / 2^1, adaptive = false, save_everystep = false)
sol2 = solve(prob, tabalg, dt = 1 / 2^1, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(prob, tabalg, dt = 1 / 2^6)
sol2 = solve(prob, BS3(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### BS5()
println("BS5")
dts = (1 / 2) .^ (6:-1:3)
sim = test_convergence(dts, probnumbig, BS5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol
sim = test_convergence(dts, probbig, BS5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol

tabalg = ExplicitRK(tableau = constructBogakiShampine5())
sol1 = solve(probnum, BS5(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnum, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(prob, BS5(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(prob, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(prob, tabalg, dt = 1 / 2^6)
sol2 = solve(prob, BS5(), dt = 1 / 2^6)

@test length(sol1) <= length(sol2) # Dual error estimators is more strict

### Tsit5()

println("Tsit5")
dts = (1 / 2) .^ (7:-1:3)
sim = test_convergence(dts, probnum, Tsit5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2
sim = test_convergence(dts, prob, Tsit5())
@test abs.(sim.𝒪est[:l2] - 5) < testTol + 0.2

tabalg = ExplicitRK(tableau = constructTsitouras5())
sol1 = solve(probnum, Tsit5(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnum, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(prob, Tsit5(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(prob, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(prob, tabalg, dt = 1 / 2^6)
sol2 = solve(prob, Tsit5(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### Vern6()

println("Vern6")
dts = (1 / 2) .^ (8:-1:5)
sim = test_convergence(dts, probnumbig, Vern6())
@test abs.(sim.𝒪est[:l2] - 6) < testTol
sim = test_convergence(dts, probbig, Vern6())
@test abs.(sim.𝒪est[:l2] - 6) < testTol

tabalg = ExplicitRK(tableau = constructVernerEfficient6(BigFloat))
sol1 = solve(probnumbig, Vern6(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(probbig, Vern6(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(probbig, tabalg, dt = 1 / 2^6)
sol2 = solve(probbig, Vern6(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### Vern7()

println("Vern7")
dts = (1 / 2) .^ (6:-1:3)
sim = test_convergence(dts, probnumbig, Vern7(), dense_errors = true)
@test abs.(sim.𝒪est[:l2] - 7) < testTol
sim = test_convergence(dts, probbig, Vern7(), dense_errors = true)
@test abs.(sim.𝒪est[:l2] - 7) < testTol

tabalg = ExplicitRK(tableau = constructVerner7(BigFloat))
sol1 = solve(probnumbig, Vern7(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(probbig, Vern7(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(probbig, tabalg, dt = 1 / 2^6)
sol2 = solve(probbig, Vern7(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### TanYam7()

println("TanYam7")
dts = (1 / 2) .^ (6:-1:3)
sim = test_convergence(dts, probnumbig, TanYam7())
@test abs.(sim.𝒪est[:l2] - 7) < testTol
sim = test_convergence(dts, probbig, TanYam7())
@test abs.(sim.𝒪est[:l2] - 7) < testTol

tabalg = ExplicitRK(tableau = constructTanakaYamashitaEfficient7(Float64))
sol1 = solve(probnum, TanYam7(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnum, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 2e-9

tabalg = ExplicitRK(tableau = constructTanakaYamashitaEfficient7(BigFloat))
sol1 = solve(probbig, TanYam7(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(prob, tabalg, dt = 1 / 2^6)
sol2 = solve(prob, TanYam7(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### Vern8()

println("Vern8")
dts = (1 / 2) .^ (6:-1:3)
sim = test_convergence(dts, probnumbig, Vern8(), dense_errors = true)
@test abs.(sim.𝒪est[:l2] - 8) < testTol
sim = test_convergence(dts, probbig, Vern8(), dense_errors = true)
@test abs.(sim.𝒪est[:l2] - 8) < testTol

tabalg = ExplicitRK(tableau = constructVerner8(BigFloat))
sol1 = solve(probnumbig, Vern8(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(probbig, Vern8(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(prob, tabalg, dt = 1 / 2^6)
sol2 = solve(prob, Vern8(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### TsitPap8()

println("TsitPap8")
dts = (1 / 2) .^ (6:-1:3)
sim = test_convergence(dts, probnumbig, TsitPap8())
@test abs.(sim.𝒪est[:l2] - 8) < testTol
sim = test_convergence(dts, probbig, TsitPap8())
@test abs.(sim.𝒪est[:l2] - 8) < testTol

tabalg = ExplicitRK(tableau = constructTsitourasPapakostas8(BigFloat))
sol1 = solve(probnumbig, TsitPap8(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test sol1.u[end] - sol2.u[end] < 1e-10

sol1 = solve(probbig, TsitPap8(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(sol1.u[end] - sol2.u[end] .< 1e-10)

sol1 = solve(prob, tabalg, dt = 1 / 2^6)
sol2 = solve(prob, TsitPap8(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)

### Vern9()

println("Vern9")
dts = (1 / 2) .^ (6:-1:3)
sim = test_convergence(dts, probnumbig, Vern9(), dense_errors = true)
@test abs.(sim.𝒪est[:l2] - 9) < testTol
sim = test_convergence(dts, probbig, Vern9(), dense_errors = true)
@test abs.(sim.𝒪est[:l2] - 9) < testTol

tabalg = ExplicitRK(tableau = constructVernerEfficient9(BigFloat))
sol1 = solve(probnumbig, Vern9(), dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg, dt = 1 / 2^6, adaptive = false, save_everystep = false)

@test abs.(sol1.u[end] - sol2.u[end]) < 1e-15

sol1 = solve(probbig, Vern9(), dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg, dt = 1 / 2^3, adaptive = false, save_everystep = false)

@test minimum(abs.(sol1.u[end] - sol2.u[end]) .< 1e-15)

sol1 = solve(probbig, tabalg, dt = 1 / 2^6)
sol2 = solve(probbig, Vern9(), dt = 1 / 2^6)

@test length(sol1) == length(sol2)
