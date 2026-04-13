using StochasticDiffEq, DiffEqNoiseProcess, Test, DiffEqDevTools, Random, LinearAlgebra
const μ = 1.01
const σ_const = 0.87

f(u, p, t) = μ * u + μ * u
f1_μ(u, p, t) = μ
f1_μ_analytic(u0, p, t, W) = u0 .* exp.((2μ - (σ_const^2) / 2)t + σ_const * W)
f2(u, p, t) = μ * u
σ(u, p, t) = σ_const * u
no_noise(u, p, t) = 0.0
f1_no_noise(u, p, t) = μ
f1_no_noise_analytic(u0, p, t, W) = u0 .* exp.(2μ * t)

ff1_μ = SplitSDEFunction(f1_μ, f2, σ, analytic = f1_μ_analytic)
prob = SDEProblem(ff1_μ, 1 / 2, (0.0, 1.0))
ff1_no_noise = SplitSDEFunction(f1_no_noise, f2, no_noise, analytic = f1_no_noise_analytic)
no_noise_prob = SDEProblem(ff1_no_noise, 1 / 2, (0.0, 1.0))

sol = solve(prob, IIF1M(), dt = 1 / 10)

prob2 = SDEProblem{false}(f, σ, 1 / 2, (0.0, 1.0), noise = NoiseWrapper(sol.W))

sol2 = solve(prob2, EM(), dt = 1 / 10)

sol = solve(no_noise_prob, IIF1M(), dt = 1 / 10)

Random.seed!(100)
dts = (1 / 2) .^ (7:-1:4) #14->7 good plot
println("IIF scalar")
sim = test_convergence(dts, prob, IIF1M(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2 # closer to 1 at this part
sim = test_convergence(dts, no_noise_prob, IIF1M(), trajectories = Int(2.0e1))
@test abs(sim.𝒪est[:l2] - 1.0) < 0.2 # closer to 1 at this part

dts = (1 / 2) .^ (7:-1:4) #14->7 good plot
println("IIF no noise scalar")
Random.seed!(100)
sim = test_convergence(dts, prob, IIF2M(), trajectories = Int(1.0e2))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2 # closer to 1 at this part
sim = test_convergence(dts, no_noise_prob, IIF2M(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 2) < 0.2 # closer to 1 at this part

#=
Random.seed!(200)
sim  = test_convergence(dts,prob,IIF1Mil(),trajectories=Int(2e1))
@test abs(sim.𝒪est[:l2]-1) < 0.3
=#

u0 = rand(2)
A = [-2.0 1.0; 1.0 -2.0]
B = [
    σ_const 0
    0 σ_const
]

function f(du, u, p, t)
    mul!(du, A, u)
    return du .+= 1.01u
end
function σ(du, u, p, t)
    mul!(@view(du[:, 1]), B, u)
    return mul!(@view(du[:, 2]), B, u)
end

function f_analytic(u0, p, t, W)
    tmp = (A + 1.01I - (B^2)) * t + B * sum(W)
    return exp(tmp) * u0
end

f1_A(du, u, p, t) = A
function f1_A_analytic(u0, p, t, W)
    tmp = (A + 1.01I - (B^2)) * t + B * sum(W)
    return exp(tmp) * u0
end
f2(du, u, p, t) = du .= μ .* u

ff1_A = SplitSDEFunction(f1_A, f2, σ, analytic = f1_A_analytic)
prob = SDEProblem(ff1_A, u0, (0.0, 1.0), noise_rate_prototype = rand(2, 2))

f1_no_noise(du, u, p, t) = A
f2(du, u, p, t) = (du .= μ .* u)
function σ22(du, u, p, t)
    return du .= 0
end
function f1_no_noise_analytic(u0, p, t, W)
    tmp = (A + 1.01I) * t
    return exp(tmp) * u0
end
ff1_A = SplitSDEFunction(f1_no_noise, f2, σ22, analytic = f1_no_noise_analytic)
prob_no_noise = SDEProblem(ff1_A, u0, (0.0, 1.0), noise_rate_prototype = rand(2, 2))

sol = solve(prob, IIF1M(), dt = 1 / 10)

dts = (1 / 2) .^ (8:-1:4) #14->7 good plot

Random.seed!(250)
println("IIF")
sim = test_convergence(dts, prob, IIF1M(), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2

sim = test_convergence(dts, prob, IIF2M(), trajectories = Int(5.0e1))
@test abs(sim.𝒪est[:l2] - 0.5) < 0.2

println("IIF no noise")
sim = test_convergence(dts, prob_no_noise, IIF1M(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 1) < 0.2

sim = test_convergence(dts, prob_no_noise, IIF2M(), trajectories = Int(1.0e1))
@test abs(sim.𝒪est[:l2] - 2) < 0.1
