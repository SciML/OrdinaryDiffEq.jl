using StochasticDiffEq, DiffEqDevTools, LinearAlgebra, Random, Test
Random.seed!(100)
dts = (1 / 2) .^ (10:-1:2) #14->7 good plot

const σ_const = 0.87
const μ = 1.01

u0 = rand(2)
A = [-2.0 1.0; 1.0 -2.0]
B = Diagonal([σ_const for i in 1:2])

function f_commute(du, u, p, t)
    mul!(du, A, u)
    return du .+= 1.01u
end

function f_commute_oop(u, p, t)
    du = A * u
    return du += 1.01u
end

function f_commute_analytic(u0, p, t, W)
    tmp = (A + 1.01I - 2 * (B^2)) * t + B * sum(W)
    return exp(tmp) * u0
end

function σ(du, u, p, t)
    du[1, 1] = σ_const * u[1]
    du[2, 1] = σ_const * u[2]
    du[1, 2] = σ_const * u[1]
    du[2, 2] = σ_const * u[2]
    du[1, 3] = σ_const * u[1]
    du[2, 3] = σ_const * u[2]
    du[1, 4] = σ_const * u[1]
    du[2, 4] = σ_const * u[2]
    return
end

function σ_oop(u, p, t)
    return σ_const * [
        u[1] u[1] u[1] u[1]
        u[2] u[2] u[2] u[2]
    ]
end

ff_commute = SDEFunction(f_commute, σ, analytic = f_commute_analytic)

prob = SDEProblem(ff_commute, u0, (0.0, 1.0), noise_rate_prototype = rand(2, 4))

sol = solve(prob, RKMilCommute(), dt = 1 / 2^(8))
sol = solve(prob, RKMilGeneral(ii_approx = IICommutative()), dt = 1 / 2^(8))
sol = solve(prob, RKMilGeneral(p = 2), dt = 1 / 2^(10))
sol = solve(prob, EM(), dt = 1 / 2^(10))

dts = (1 / 2) .^ (10:-1:3) #14->7 good plot
sim2 = test_convergence(dts, prob, EM(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:final] - 0.5) < 0.2

sim2 = test_convergence(dts, prob, RKMilCommute(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:final] - 1.0) < 0.2

sim2 = test_convergence(dts, prob, RKMilGeneral(ii_approx = IICommutative()), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:final] - 1.0) < 0.2

sim2 = test_convergence(dts, prob, RKMilGeneral(p = 2), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:final] - 1.0) < 0.2

ff_commute_oop = SDEFunction(f_commute_oop, σ_oop, analytic = f_commute_analytic)

proboop = SDEProblem(ff_commute_oop, u0, (0.0, 1.0), noise_rate_prototype = rand(2, 4))

sol = solve(proboop, RKMilCommute(), dt = 1 / 2^(8))
sol = solve(proboop, RKMilGeneral(ii_approx = IICommutative()), dt = 1 / 2^(8))
sol = solve(proboop, RKMilGeneral(p = 2), dt = 1 / 2^(10))
sol = solve(proboop, EM(), dt = 1 / 2^(10))

dts = (1 / 2) .^ (10:-1:3) #14->7 good plot
sim2 = test_convergence(dts, proboop, EM(), trajectories = Int(1.0e3))
@test abs(sim2.𝒪est[:final] - 0.5) < 0.2

sim2 = test_convergence(dts, proboop, RKMilCommute(), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:final] - 1.0) < 0.2

sim2 = test_convergence(dts, proboop, RKMilGeneral(ii_approx = IICommutative()), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:final] - 1.0) < 0.2

sim2 = test_convergence(dts, proboop, RKMilGeneral(p = 2), trajectories = Int(2.0e2))
@test abs(sim2.𝒪est[:final] - 1.0) < 0.2
