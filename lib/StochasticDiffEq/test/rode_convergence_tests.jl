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

const TAYLOR_STEP_COUNTS = STEP_COUNTS[1:5]

function taylor15_study(alg; seed = 20260921, sigma = 1.0)
    rng = MersenneTwister(seed)
    errors = zeros(PATHS, length(TAYLOR_STEP_COUNTS))
    floors = zeros(PATHS)
    for m in 1:PATHS
        path = sigma .* wiener_path(rng)
        g = cos.(5 .* path)
        exact = exp.(-cumulative_trapezoid(g, FINE_DT))
        floors[m] = maximum(abs, exp.(-cumulative_trapezoid(g[1:4:end], 4FINE_DT)) .- exact[1:4:end])
        noise = NoiseGrid(FINE_GRID, path)
        for (j, n) in enumerate(TAYLOR_STEP_COUNTS)
            sol = decay_solution(alg, noise, n, 1.0)
            stride = FINE_POINTS ÷ n
            errors[m, j] = maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n)
        end
    end
    strong = [sqrt(mean(errors[:, j] .^ 2)) for j in eachindex(TAYLOR_STEP_COUNTS)]
    return (
        strong = strong, strong_order = fitted_slope(TAYLOR_STEP_COUNTS, strong),
        reference_floor = maximum(floors),
    )
end

@testset "RandomTaylor15 reaches order 1.5 on a resolved path" begin
    taylor = taylor15_study(RandomTaylor15())
    euler = taylor15_study(RandomEM())
    @test minimum(taylor.strong) > 5 * taylor.reference_floor
    @test 1.5 < taylor.strong_order < 2.2
    @test taylor.strong_order > euler.strong_order + 0.5
    @test taylor.strong[end] < euler.strong[end] / 10
end

@testset "RandomTaylor15 keeps its order on a Brownian path of another amplitude" begin
    taylor = taylor15_study(RandomTaylor15(); sigma = 0.2)
    euler = taylor15_study(RandomEM(); sigma = 0.2)
    @test minimum(taylor.strong) > 3 * taylor.reference_floor
    @test 1.5 < taylor.strong_order < 2.3
    @test taylor.strong[end] < euler.strong[end] / 20
end

@testset "RandomTaylor15 is Heun in the deterministic limit" begin
    smooth(u, p, t, W) = -u
    exact = exp.(-FINE_GRID)
    noise = NoiseGrid(FINE_GRID, zeros(FINE_POINTS + 1))
    function smooth_solution(alg, n)
        prob = RODEProblem{false}(smooth, 1.0, (0.0, TEND), noise = noise)
        return solve(prob, alg, dt = TEND / n, save_everystep = true, adaptive = false)
    end
    function smooth_error(alg, n)
        sol = smooth_solution(alg, n)
        stride = FINE_POINTS ÷ n
        return maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n)
    end
    errors = [smooth_error(RandomTaylor15(), n) for n in STEP_COUNTS]
    @test abs(fitted_slope(STEP_COUNTS, errors) - 2) < 0.2
    taylor = smooth_solution(RandomTaylor15(), STEP_COUNTS[end]).u
    heun = smooth_solution(RandomHeun(), STEP_COUNTS[end]).u
    @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(taylor, heun))
end

@testset "RandomTaylor15 in-place matches out-of-place" begin
    noise = NoiseGrid(FINE_GRID, wiener_path(MersenneTwister(20260921)))
    u0 = [1.0, 2.0]
    for n in (TAYLOR_STEP_COUNTS[1], TAYLOR_STEP_COUNTS[end])
        outofplace = decay_solution(RandomTaylor15(), noise, n, u0).u
        inplace = decay_solution(RandomTaylor15(), noise, n, copy(u0); inplace = true).u
        @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(outofplace, inplace))
    end
end

@testset "RandomTaylor15 step integrals are exact on a misaligned grid" begin
    grid = collect(range(0.0, 1.0; length = 13))
    integrals = StochasticDiffEq.StochasticDiffEqRODE.path_integrals
    for (t, dt) in ((0.0, 0.25), (0.05, 0.25), (0.03, 0.04), (1 / 12, 1 / 12), (0.5, 0.5))
        W = (t = grid, W = copy(grid), dW = dt, curW = t)
        I1, I2 = integrals(W, t, dt, t)
        @test I1 ≈ dt^2 / 2
        @test I2 ≈ dt^3 / 3
    end
end

@testset "RandomTaylor15 integrates backwards in time" begin
    n = 64
    back_grid = collect(range(TEND, 0.0; length = FINE_POINTS + 1))
    path = reverse(wiener_path(MersenneTwister(20260921)))
    g = cos.(5 .* path)
    exact = exp.(-cumulative_trapezoid(g, -FINE_DT))
    prob = RODEProblem{false}(decay_oop, 1.0, (TEND, 0.0), noise = NoiseGrid(back_grid, path))
    sol = solve(prob, RandomTaylor15(), dt = -TEND / n, save_everystep = true, adaptive = false)
    stride = FINE_POINTS ÷ n
    @test maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n) < 1.0e-2
end

@testset "RandomTaylor15 rejects noise it cannot resolve" begin
    prob = RODEProblem{false}(decay_oop, 1.0, (0.0, TEND))
    @test_throws "not compatible with the chosen noise type" solve(
        prob, RandomTaylor15(), dt = TEND / 16, adaptive = false
    )
end

using LinearAlgebra

const AFFINE_POINTS = 2^13
const AFFINE_GRID = collect(range(0.0, TEND; length = AFFINE_POINTS + 1))
const AFFINE_STEP_COUNTS = [16, 32, 64, 128, 256]
const AFFINE_PATHS = 8
const AFFINE_U0 = [0.7, -0.4]

affine_scalar(u, p, t, W) = -u + sin(3t) + (cos(u) + t / 2) * W[1]
function affine_scalar!(du, u, p, t, W)
    du[1] = affine_scalar(u[1], p, t, W)
    return nothing
end

function affine_pair(u, p, t, W)
    return [
        -u[1] + sin(u[2]) + t + u[2] * W[1] + (1 + t) * cos(u[1]) * W[2],
        cos(2t) * u[1] - u[2]^2 / 2 + (0.3 * u[1] * u[2] + t) * W[1] +
            (u[1]^2 - u[2]) * W[2],
    ]
end
function affine_pair!(du, u, p, t, W)
    du .= affine_pair(u, p, t, W)
    return nothing
end

function affine_reference(f, u0, grid, path; substeps = 4)
    u = u0
    reference = [u0]
    for k in 1:(length(grid) - 1)
        h = (grid[k + 1] - grid[k]) / substeps
        dw = (path[k + 1] - path[k]) / substeps
        for s in 0:(substeps - 1)
            t = grid[k] + s * h
            w = path[k] + s * dw
            k1 = f(u, nothing, t, w)
            k2 = f(u + (h / 2) * k1, nothing, t + h / 2, w + dw / 2)
            k3 = f(u + (h / 2) * k2, nothing, t + h / 2, w + dw / 2)
            k4 = f(u + h * k3, nothing, t + h, w + dw)
            u = u + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        end
        push!(reference, u)
    end
    return reference
end

function affine_solution(f, u0, grid, path, alg, n)
    tspan = (grid[1], grid[end])
    prob = RODEProblem(f, u0, tspan, noise = NoiseGrid(grid, path))
    dt = (tspan[2] - tspan[1]) / n
    return solve(prob, alg, dt = dt, save_everystep = true, adaptive = false)
end

function affine_errors(f, u0, grid, path, alg, step_counts)
    reference = affine_reference(f, u0, grid, path)
    return map(step_counts) do n
        sol = affine_solution(f, u0, grid, path, alg, n)
        stride = (length(grid) - 1) ÷ n
        maximum(maximum(abs, sol.u[k + 1] - reference[k * stride + 1]) for k in 0:n)
    end
end

function affine_study(f, u0, grid, paths, alg, step_counts)
    errors = [affine_errors(f, u0, grid, path, alg, step_counts) for path in paths]
    strong = [sqrt(mean(e[j]^2 for e in errors)) for j in eachindex(step_counts)]
    return (strong = strong, order = fitted_slope(step_counts, strong))
end

function affine_wiener(rng, offset)
    W = zeros(AFFINE_POINTS + 1)
    W[2:end] .= cumsum(sqrt(TEND / AFFINE_POINTS) .* randn(rng, AFFINE_POINTS))
    return W .+ offset
end

function affine_wiener_pair(rng)
    W1 = affine_wiener(rng, 0.5)
    W2 = affine_wiener(rng, -0.3)
    return [[W1[k], W2[k]] for k in eachindex(W1)]
end

function fractional_paths(rng, hurst, points, paths)
    t = (1:points) .* (TEND / points)
    covariance = [(s^(2hurst) + r^(2hurst) - abs(s - r)^(2hurst)) / 2 for s in t, r in t]
    L = cholesky(Symmetric(covariance)).L
    return [[0.0; L * randn(rng, points)] for _ in 1:paths]
end

@testset "RandomAffineHeun reaches order 2 on a Wiener path" begin
    rng = MersenneTwister(20261004)
    paths = [affine_wiener(rng, 0.0) for _ in 1:AFFINE_PATHS]
    affine = affine_study(
        affine_scalar, 0.5, AFFINE_GRID, paths, RandomAffineHeun(), AFFINE_STEP_COUNTS
    )
    @test abs(affine.order - 2) < 0.2
end

@testset "RandomAffineHeun reaches order 2 with two noises that do not commute" begin
    rng = MersenneTwister(20261004)
    paths = [affine_wiener_pair(rng) for _ in 1:AFFINE_PATHS]
    affine = affine_study(
        affine_pair, AFFINE_U0, AFFINE_GRID, paths, RandomAffineHeun(),
        2 .* AFFINE_STEP_COUNTS
    )
    @test abs(affine.order - 2) < 0.2
end

@testset "RandomAffineHeun keeps the two orderings of a double integral apart" begin
    path = [[sin(128pi * t), cos(128pi * t)] for t in AFFINE_GRID]
    errors = affine_errors(
        affine_pair, AFFINE_U0, AFFINE_GRID, path, RandomAffineHeun(), AFFINE_STEP_COUNTS
    )
    @test abs(fitted_slope(AFFINE_STEP_COUNTS, errors) - 2) < 0.3
    @test errors[3] < 3.0e-4
end

@testset "RandomAffineHeun keeps order 2 on a fractional Brownian path, Hurst index 0.3" begin
    points = 2^11
    grid = collect(range(0.0, TEND; length = points + 1))
    step_counts = [8, 16, 32, 64, 128]
    paths = fractional_paths(MersenneTwister(20261004), 0.3, points, AFFINE_PATHS)
    affine = affine_study(affine_scalar, 0.5, grid, paths, RandomAffineHeun(), step_counts)
    heun = affine_study(affine_scalar, 0.5, grid, paths, RandomHeun(), step_counts)
    euler = affine_study(affine_scalar, 0.5, grid, paths, RandomEM(), step_counts)
    @test abs(affine.order - 2) < 0.2
    @test heun.order < 1.3
    @test euler.order < 1.3
    @test affine.strong[end] < heun.strong[end] / 50
end

@testset "RandomAffineHeun is Heun when the noise leaves the right-hand side" begin
    drift(u, p, t, W) = -u + sin(3t)
    path = affine_wiener(MersenneTwister(20261004), 0.5)
    exact = @. 0.8 * exp(-AFFINE_GRID) + (sin(3 * AFFINE_GRID) - 3 * cos(3 * AFFINE_GRID)) / 10
    errors = map(AFFINE_STEP_COUNTS) do n
        sol = affine_solution(drift, 0.5, AFFINE_GRID, path, RandomAffineHeun(), n)
        stride = AFFINE_POINTS ÷ n
        maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n)
    end
    @test abs(fitted_slope(AFFINE_STEP_COUNTS, errors) - 2) < 0.2
    affine = affine_solution(drift, 0.5, AFFINE_GRID, path, RandomAffineHeun(), 64).u
    heun = affine_solution(drift, 0.5, AFFINE_GRID, path, RandomHeun(), 64).u
    @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(affine, heun))
end

@testset "RandomAffineHeun in-place matches out-of-place" begin
    pair = affine_wiener_pair(MersenneTwister(20261004))
    alg = RandomAffineHeun()
    for (grid, path) in ((AFFINE_GRID, pair), (reverse(AFFINE_GRID), reverse(pair)))
        single = first.(path)
        for n in (16, 100, 256)
            outofplace = affine_solution(affine_pair, AFFINE_U0, grid, path, alg, n).u
            inplace = affine_solution(affine_pair!, copy(AFFINE_U0), grid, path, alg, n).u
            @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(outofplace, inplace))
            scalar = affine_solution(affine_scalar, 0.5, grid, single, alg, n).u
            inplace = affine_solution(affine_scalar!, [0.5], grid, single, alg, n).u
            vector = affine_solution(affine_scalar, 0.5, grid, [[w] for w in single], alg, n).u
            @test all(isapprox(a, b[1], rtol = 1.0e-9) for (a, b) in zip(scalar, inplace))
            @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(scalar, vector))
        end
    end
end

@testset "RandomAffineHeun integrates backwards in time" begin
    rng = MersenneTwister(20261004)
    paths = [reverse(affine_wiener_pair(rng)) for _ in 1:AFFINE_PATHS]
    affine = affine_study(
        affine_pair, AFFINE_U0, reverse(AFFINE_GRID), paths, RandomAffineHeun(),
        2 .* AFFINE_STEP_COUNTS
    )
    @test abs(affine.order - 2) < 0.2
end

@testset "RandomAffineHeun converges with steps finer than the grid" begin
    rng = MersenneTwister(20261004)
    grid = collect(range(0.0, TEND; length = 101))
    path = [[0.5 + randn(rng) / 2, -0.3 + randn(rng) / 2] for _ in grid]
    step_counts = [384, 768, 1536]
    for (grid, path) in ((grid, path), (reverse(grid), reverse(path)))
        reference = affine_reference(affine_pair, AFFINE_U0, grid, path; substeps = 256)[end]
        errors = map(step_counts) do n
            sol = affine_solution(affine_pair!, copy(AFFINE_U0), grid, path, RandomAffineHeun(), n)
            maximum(abs, sol.u[end] - reference)
        end
        @test abs(fitted_slope(step_counts, errors) - 2) < 0.3
        @test errors[end] < 1.0e-6
    end
end

@testset "RandomAffineHeun step integrals are exact on a misaligned grid" begin
    integrals! = StochasticDiffEq.StochasticDiffEqRODE.affine_integrals!
    value(t) = t < 1 / 2 ? [1 / 2 - t, 1 + 2t] : [t - 1 / 2, 7 / 2 - 3t]
    slope(t) = t < 1 / 2 ? [-1.0, 2.0] : [1.0, -3.0]
    function piece(t, h)
        a, b = value(t), slope(t + h / 2)
        I1 = @. a * h + b * h^2 / 2
        I01 = @. a * h^2 / 2 + b * h^3 / 3
        I2 = [
            a[i] * a[j] * h^2 / 2 + b[i] * a[j] * h^3 / 6 + a[i] * b[j] * h^3 / 3 +
                b[i] * b[j] * h^4 / 8 for i in 1:2, j in 1:2
        ]
        return I1, I01, I2
    end
    ascending = collect(range(0.0, 1.0; length = 13))
    forward = (
        (0.0, 0.25), (0.05, 0.25), (0.03, 0.04), (0.3, 0.45), (5 / 12, 1 / 6), (0.45, 0.1),
        (0.5, 0.5),
    )
    backward = (
        (1.0, -0.25), (0.8, -0.45), (0.52, -0.04), (0.97, -0.04), (7 / 12, -1 / 6),
        (0.5, -0.5),
    )
    for (grid, steps) in ((ascending, forward), (reverse(ascending), backward))
        for (t, dt) in steps
            kink = clamp(1 / 2, minmax(t, t + dt)...)
            A1, A01, A2 = piece(t, kink - t)
            B1, B01, B2 = piece(kink, t + dt - kink)
            W = (t = grid, W = value.(grid), curW = value(t))
            I1, I01, I2 = zeros(2), zeros(2), zeros(2, 2)
            integrals!(I1, I01, I2, W, t, dt, value(t + dt))
            @test I1 ≈ A1 + B1 atol = 1.0e-14
            @test I01 ≈ A01 + B01 + (kink - t) * B1 atol = 1.0e-14
            @test I2 ≈ A2 + B2 + A1 * B1' atol = 1.0e-14
        end
    end
end

@testset "RandomAffineHeun rejects noise that is not a NoiseGrid of scalars or vectors" begin
    prob = RODEProblem(affine_scalar, 0.5, (0.0, TEND))
    @test_throws "not compatible with the chosen noise type" solve(
        prob, RandomAffineHeun(), dt = TEND / 16, adaptive = false
    )
    noise = NoiseGrid(AFFINE_GRID, [fill(t, 1, 1) for t in AFFINE_GRID])
    prob = RODEProblem(affine_scalar, 0.5, (0.0, TEND), noise = noise)
    @test_throws "not compatible with the chosen noise type" solve(
        prob, RandomAffineHeun(), dt = TEND / 16, adaptive = false
    )
end
