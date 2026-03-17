using StochasticDiffEq, DiffEqDevTools, LinearAlgebra, Random, Test
Random.seed!(100)
dts = (1 / 2) .^ (10:-1:2) #14->7 good plot

# From RUNGE–KUTTA METHODS FOR THE STRONG APPROXIMATION OF SOLUTIONS OF STOCHASTIC DIFFERENTIAL EQUATIONS
# (7.4)

d = 4;
m = 10
u0 = rand(4)
A = zeros(d, d)
for i in 1:d, j in 1:d

    global A
    i == j && (A[i, j] = -3 / 2)
    i != j && (A[i, j] = 1 / 20)
end

B = [zeros(d, d) for i in 1:m]
for k in 1:m, i in 1:d, j in 1:d
    global B
    i == j && (B[k][i, j] = 1 / 5)
    i != j && (B[k][i, j] = 1 / 100)
end

function f_noncommute(du, u, p, t)
    return mul!(du, A, u)
end

function g_noncommute(du, u, p, t)
    for i in 1:m
        mul!(@view(du[:, i]), B[i], u)
    end
    return
end

function f_noncommute_analytic(u0, p, t, W)
    tmp = (A - 0.5 * sum(B[i]^2 for i in 1:m)) * t + sum(B[i] * W[i] for i in 1:m)
    return exp(tmp) * u0
end

function f_noncommute_analytic_stratonovich(u0, p, t, W)
    tmp = A * t + sum(B[i] * W[i] for i in 1:m)
    return exp(tmp) * u0
end

ff_noncommute = SDEFunction(f_noncommute, g_noncommute, analytic = f_noncommute_analytic)

prob = SDEProblem(ff_noncommute, u0, (0.0, 1.0), noise_rate_prototype = rand(4, m))

sol = solve(prob, EM(), dt = 1 / 2^(8))
sol = solve(prob, RKMilGeneral(p = 10), dt = 1 / 2^(8))

dts = (1 / 2) .^ (10:-1:3) #14->7 good plot
sim1 = test_convergence(dts, prob, EM(), trajectories = Int(1.0e2))
@test abs(sim1.𝒪est[:final] - 0.5) < 0.2
sim2 = test_convergence(dts, prob, RKMilCommute(), trajectories = Int(1.0e2))
@test abs(sim2.𝒪est[:final] - 1) < 0.2
sim3 = test_convergence(dts, prob, RKMilGeneral(p = 2), trajectories = Int(1.0e2))
@test abs(sim3.𝒪est[:final] - 1) < 0.2

ff_noncommute_stratonovich = SDEFunction(
    f_noncommute, g_noncommute, analytic = f_noncommute_analytic_stratonovich
)
prob_stratonovich = SDEProblem(ff_noncommute_stratonovich, u0, (0.0, 1.0), noise_rate_prototype = rand(4, m))

sim4 = test_convergence(dts, prob_stratonovich, EulerHeun(), trajectories = Int(1.0e2))
@test abs(sim4.𝒪est[:final] - 1.0) < 0.2

d = 2;
m = 4
u0 = [2.0, 2.0]
α = 1 / 2

function f_noncommute_2(du, u, p, t)
    du .= 0
    return nothing
end

function g_noncommute_2(du, u, p, t)
    du[1, 1] = cos(p[1]) * sin(u[1])
    du[2, 1] = sin(p[1]) * sin(u[1])
    du[1, 2] = cos(p[1]) * cos(u[1])
    du[2, 2] = sin(p[1]) * cos(u[1])
    du[1, 3] = -sin(p[1]) * sin(u[2])
    du[2, 3] = cos(p[1]) * sin(u[2])
    du[1, 4] = -sin(p[1]) * cos(u[2])
    du[2, 4] = cos(p[1]) * cos(u[2])
    return nothing
end

p = [α]
prob2 = SDEProblem(
    f_noncommute_2, g_noncommute_2, u0, (0.0, 0.5), p, noise_rate_prototype = rand(2, m)
)

sol1 = solve(prob2, EM(), dt = 1 / 2^(8))
sol2 = solve(prob2, RKMilGeneral(p = true, dt = 1 / 2^(7)), dt = 1 / 2^(7), adaptive = false)
sol3 = solve(prob2, RKMilGeneral(), dt = 1 / 2^(7))

dts = (1 / 2) .^ (6:-1:2) #14->7 good plot
test_dt = 1 / 2^(10)
sim5 = analyticless_test_convergence(
    dts, prob2, EM(), test_dt, trajectories = 300, use_noise_grid = false
)
@test abs(sim5.𝒪est[:final] - 0.5) < 0.2
sim6 = analyticless_test_convergence(
    dts, prob2, RKMilGeneral(p = true, dt = test_dt),
    test_dt, trajectories = 100, use_noise_grid = false
)
@test_broken abs(sim6.𝒪est[:final] - 1.0) < 0.2
@test abs(sim6.𝒪est[:weak_final] - 1.0) < 0.2
sim7 = analyticless_test_convergence(
    dts, prob2, EulerHeun(), test_dt, trajectories = 300, use_noise_grid = false
)
@test abs(sim7.𝒪est[:final] - 0.5) < 0.2
