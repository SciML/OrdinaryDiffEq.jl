using Random, StochasticDiffEq, DiffEqDevTools, Test

# dX = μ dt + s dW has X(T) ~ Normal(X0 + μT, s²T), so the 95% half-width of a
# sample mean over N paths is 1.96 s √T / √N.
Random.seed!(100)
s, T = [0.5, 2.0], 1.0
brownian = SDEFunction(
    (u, p, t) -> fill!(similar(u), 1.0), (u, p, t) -> s;
    analytic = (u0, p, t, W) -> u0 .+ t .+ s .* W
)
Ns = [10, 100, 500]
setup = Dict(:alg => EM(), :dts => [T / 4], :adaptive => false)

scalar_prob = SDEProblem(
    SDEFunction(
        (u, p, t) -> 1.0, (u, p, t) -> s[1];
        analytic = (u0, p, t, W) -> u0 + t + s[1] * W
    ), 0.0, (0.0, T)
)
halfwidth = 1.96 * s[1] * sqrt(T) ./ sqrt.(Ns)
se = get_sample_errors(scalar_prob, setup, numruns = Ns, solution_runs = 200)
@test all(0.8 .< se ./ halfwidth .< 1.2)
se100 = get_sample_errors(scalar_prob, setup, numruns = 100, solution_runs = 200)
@test 0.8 < se100 / halfwidth[2] < 1.2

vector_prob = SDEProblem(brownian, zeros(2), (0.0, T))
halfwidth_vec = 1.96 * sqrt(sum(abs2, s) * T) ./ sqrt.(Ns)
se_vec = get_sample_errors(vector_prob, setup, numruns = Ns, solution_runs = 200)
@test all(0.8 .< se_vec ./ halfwidth_vec .< 1.2)
