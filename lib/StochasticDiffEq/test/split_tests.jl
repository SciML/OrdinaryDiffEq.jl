using StochasticDiffEq, DiffEqNoiseProcess, Test, DiffEqDevTools

f(u, p, t) = (1.01) * u
f1(u, p, t) = (1.01) / 2 * u
f2(u, p, t) = (1.01) / 2 * u
σ(u, p, t) = 0.87u
f_split_analytic(u0, p, t, W) = @. u0 * exp(0.63155t + 0.87W)

prob = SplitSDEProblem{false}(f1, f2, σ, 1 / 2, (0.0, 1.0))
sol = solve(prob, SplitEM(), dt = 1 / 10, save_noise = true)

prob = SDEProblem{false}(f, σ, 1 / 2, (0.0, 1.0), noise = NoiseWrapper(sol.W))
sol2 = solve(prob, EM(), dt = 1 / 10)

@test sol.u ≈ sol2.u

u0 = rand(4)

ff_split = SplitSDEFunction(f1, f2, σ, analytic = f_split_analytic)
prob = SplitSDEProblem(ff_split, u0, (0.0, 1.0))

sol = solve(prob, SplitEM(), dt = 1 / 10, save_noise = true)

prob = SDEProblem{false}(f, σ, u0, (0.0, 1.0), noise = NoiseWrapper(sol.W))

sol2 = solve(prob, EM(), dt = 1 / 10)

@test sol.u[end][:] ≈ sol2.u[end][:]

################################################################################

### Only first

println("Only First")

α = 0.1
β = 0.5
ff1 = (u, p, t) -> β ./ sqrt.(1 + t) - u ./ (2 * (1 + t))
ff2 = (u, p, t) -> 0.0
σ2 = (u, p, t) -> α * β ./ sqrt.(1 + t)
ff1_analytic(u0, p, t, W) = @. u0 / sqrt(1 + t) + β * (t + α * W) / sqrt(1 + t)
f_ff1 = SplitSDEFunction(ff1, ff2, σ2, analytic = ff1_analytic)
prob = SplitSDEProblem(f_ff1, 1.0, (0.0, 1.0))

sol = solve(prob, EM(), dt = 1 / 10)
sol2 = solve(prob, SKenCarp(), dt = 1 / 10)

dts = (1 / 2) .^ (10:-1:2) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3

### Only second

println("Only Second")

α = 0.1
β = 0.5
ff1 = (u, p, t) -> 0.0
ff2 = (u, p, t) -> β ./ sqrt.(1 + t) - u ./ (2 * (1 + t))
σ2 = (u, p, t) -> α * β ./ sqrt.(1 + t)
f_ff1 = SplitSDEFunction(ff1, ff2, σ2, analytic = ff1_analytic)
prob = SplitSDEProblem(f_ff1, 1.0, (0.0, 1.0))

sol = solve(prob, EM(), dt = 1 / 10)
sol2 = solve(prob, SKenCarp(), dt = 1 / 10, seed = 1)

dts = (1 / 2) .^ (10:-1:2) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3

### Both

println("Both")

α = 0.1
β = 0.5
ff1 = (u, p, t) -> β ./ sqrt.(1 + t)
ff2 = (u, p, t) -> - u ./ (2 * (1 + t))
σ2 = (u, p, t) -> α * β ./ sqrt.(1 + t)
f_ff1 = SplitSDEFunction(ff1, ff2, σ2, analytic = ff1_analytic)
prob = SplitSDEProblem(f_ff1, 1.0, (0.0, 1.0))

sol = solve(prob, EM(), dt = 1 / 10)
sol2 = solve(prob, SKenCarp(), dt = 1 / 10)

dts = (1 / 2) .^ (10:-1:2) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3
sim10 = test_convergence(
    dts, prob, SKenCarp(nlsolve = StochasticDiffEq.NLFunctional()), trajectories = Int(1.0e1)
)
@test abs(sim10.𝒪est[:final] - 2) < 0.3

################################################################################

### Only first

println("Only First")

α = 0.1
β = 0.5
ff1 = (du, u, p, t) -> @. du = β / sqrt(1 + t) - u / (2 * (1 + t))
ff2 = (du, u, p, t) -> @. du = 0.0
σ2 = (du, u, p, t) -> @. du = α * β / sqrt(1 + t)
f_ff1 = SplitSDEFunction(ff1, ff2, σ2, analytic = ff1_analytic)
prob = SplitSDEProblem(f_ff1, [1.0], (0.0, 1.0))

sol = solve(prob, EM(), dt = 1 / 10)
sol2 = solve(prob, SKenCarp(), dt = 1 / 10)

dts = (1 / 2) .^ (10:-1:2) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3

### Only second

println("Only Second")

α = 0.1
β = 0.5
ff1 = (du, u, p, t) -> @. du = 0.0
ff2 = (du, u, p, t) -> @. du = β / sqrt(1 + t) - u / (2 * (1 + t))
σ2 = (du, u, p, t) -> @. du = α * β / sqrt(1 + t)
f_ff1 = SplitSDEFunction(ff1, ff2, σ2, analytic = ff1_analytic)
prob = SplitSDEProblem(f_ff1, [1.0], (0.0, 1.0))

sol = solve(prob, EM(), dt = 1 / 10)
sol2 = solve(prob, SKenCarp(), dt = 1 / 10)

dts = (1 / 2) .^ (10:-1:2) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3

### Both

println("Both")

α = 0.1
β = 0.5
ff1 = (du, u, p, t) -> @. du = β / sqrt(1 + t)
ff2 = (du, u, p, t) -> @. du = - u / (2 * (1 + t))
σ2 = (du, u, p, t) -> @. du = α * β / sqrt(1 + t)
f_ff1 = SplitSDEFunction(ff1, ff2, σ2, analytic = ff1_analytic)
prob = SplitSDEProblem(f_ff1, [1.0], (0.0, 1.0))

sol = solve(prob, EM(), dt = 1 / 10)
sol2 = solve(prob, SKenCarp(), dt = 1 / 10)

dts = (1 / 2) .^ (10:-1:2) #14->7 good plot
sim10 = test_convergence(dts, prob, SKenCarp(), trajectories = Int(1.0e1))
@test abs(sim10.𝒪est[:final] - 2) < 0.3
sim10 = test_convergence(
    dts, prob, SKenCarp(nlsolve = StochasticDiffEq.NLFunctional()), trajectories = Int(1.0e1)
)
@test abs(sim10.𝒪est[:final] - 2) < 0.3
