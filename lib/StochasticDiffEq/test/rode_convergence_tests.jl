using StochasticDiffEq, DiffEqNoiseProcess, Random, Statistics, Test

const TEND = 1.0
const FINE_POINTS = 2^19
const FINE_DT = TEND / FINE_POINTS
const FINE_GRID = collect(range(0.0, TEND; length = FINE_POINTS + 1))
const STEP_COUNTS = [FINE_POINTS ÷ 2^k for k in 15:-1:9]
const PATHS = 32
const ALGORITHMS = (RandomEM(), RandomHeun(), RandomTamedEM())

function wiener_path(rng)
    W = zeros(FINE_POINTS + 1)
    W[2:end] .= cumsum(sqrt(FINE_DT) .* randn(rng, FINE_POINTS))
    return W
end

function cumulative_trapezoid(g, dt)
    I = zeros(length(g))
    acc = zero(eltype(g))
    @inbounds for k in 1:(length(g) - 1)
        acc += dt * (g[k] + g[k + 1]) / 2
        I[k + 1] = acc
    end
    return I
end

function fitted_slope(step_counts, errors)
    x = log2.(1 ./ step_counts)
    y = log2.(errors)
    n = length(x)
    return (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x .^ 2) - sum(x)^2)
end

decay_oop(u, p, t, W) = -u .* cos(5W)
decay_iip(du, u, p, t, W) = (du .= -u .* cos(5W))

function decay_solution(alg, noise, n, u0; inplace = false)
    prob = if inplace
        RODEProblem(decay_iip, u0, (0.0, TEND), noise = noise)
    else
        RODEProblem{false}(decay_oop, u0, (0.0, TEND), noise = noise)
    end
    return solve(prob, alg, dt = TEND / n, save_everystep = true, adaptive = false)
end

function decay_study(alg; seed = 20260918)
    rng = MersenneTwister(seed)
    errors = zeros(PATHS, length(STEP_COUNTS))
    floors = zeros(PATHS)
    for m in 1:PATHS
        path = wiener_path(rng)
        g = cos.(5 .* path)
        exact = exp.(-cumulative_trapezoid(g, FINE_DT))
        floors[m] = maximum(abs, exp.(-cumulative_trapezoid(g[1:4:end], 4FINE_DT)) .- exact[1:4:end])
        noise = NoiseGrid(FINE_GRID, path)
        for (j, n) in enumerate(STEP_COUNTS)
            sol = decay_solution(alg, noise, n, 1.0)
            stride = FINE_POINTS ÷ n
            errors[m, j] = maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n)
        end
    end
    strong = [sqrt(mean(errors[:, j] .^ 2)) for j in eachindex(STEP_COUNTS)]
    pathwise = [fitted_slope(STEP_COUNTS, errors[m, :]) for m in 1:PATHS]
    return (
        strong = strong, strong_order = fitted_slope(STEP_COUNTS, strong),
        pathwise_order = median(pathwise), reference_floor = maximum(floors),
    )
end

@testset "classical orders are visible when the noise leaves the right-hand side" begin
    smooth(u, p, t, W) = -u
    exact = exp.(-FINE_GRID)
    noise = NoiseGrid(FINE_GRID, zeros(FINE_POINTS + 1))
    function smooth_error(alg, n)
        prob = RODEProblem{false}(smooth, 1.0, (0.0, TEND), noise = noise)
        sol = solve(prob, alg, dt = TEND / n, save_everystep = true, adaptive = false)
        stride = FINE_POINTS ÷ n
        return maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n)
    end
    @test abs(fitted_slope(STEP_COUNTS, [smooth_error(RandomEM(), n) for n in STEP_COUNTS]) - 1) < 0.2
    @test abs(fitted_slope(STEP_COUNTS, [smooth_error(RandomHeun(), n) for n in STEP_COUNTS]) - 2) < 0.2
end

@testset "Wiener-driven RODE order: $(nameof(typeof(alg)))" for alg in ALGORITHMS
    study = decay_study(alg)
    @test minimum(study.strong) > 20 * study.reference_floor
    @test abs(study.strong_order - 1) < 0.2
    @test abs(study.pathwise_order - 1) < 0.2
end

@testset "in-place matches out-of-place: $(nameof(typeof(alg)))" for alg in ALGORITHMS
    noise = NoiseGrid(FINE_GRID, wiener_path(MersenneTwister(20260918)))
    u0 = [1.0, 1.0]
    for n in (STEP_COUNTS[1], STEP_COUNTS[4], STEP_COUNTS[end])
        @test decay_solution(alg, noise, n, u0).u ==
            decay_solution(alg, noise, n, copy(u0); inplace = true).u
    end
end
