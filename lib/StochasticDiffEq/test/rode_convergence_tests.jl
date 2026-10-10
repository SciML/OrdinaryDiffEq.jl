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

function cosine_decay_exact(V)
    exact = ones(length(V))
    acc = 0.0
    for k in 1:(length(V) - 1)
        a = V[k]
        b = V[k + 1]
        acc += abs(b - a) < 1.0e-8 ? FINE_DT * cos((a + b) / 2) :
            FINE_DT * (sin(b) - sin(a)) / (b - a)
        exact[k + 1] = exp(-acc)
    end
    return exact
end

function vector_decay_study(alg, c; seed = 20261004)
    rng = MersenneTwister(seed)
    m = length(c)
    decay(u, p, t, W) = -u * cos(sum(c .* W))
    errors = zeros(PATHS, length(TAYLOR_STEP_COUNTS))
    for q in 1:PATHS
        paths = [wiener_path(rng) for _ in 1:m]
        exact = cosine_decay_exact(sum(c[i] .* paths[i] for i in 1:m))
        noise = NoiseGrid(FINE_GRID, [[paths[i][k] for i in 1:m] for k in 1:(FINE_POINTS + 1)])
        for (j, n) in enumerate(TAYLOR_STEP_COUNTS)
            prob = RODEProblem{false}(decay, 1.0, (0.0, TEND), noise = noise)
            sol = solve(prob, alg, dt = TEND / n, save_everystep = true, adaptive = false)
            stride = FINE_POINTS ÷ n
            errors[q, j] = maximum(abs(sol.u[k + 1] - exact[k * stride + 1]) for k in 0:n)
        end
    end
    strong = [sqrt(mean(errors[:, j] .^ 2)) for j in eachindex(TAYLOR_STEP_COUNTS)]
    return (strong = strong, strong_order = fitted_slope(TAYLOR_STEP_COUNTS, strong))
end

@testset "RandomTaylor15 keeps its order with $(length(c)) noise components" for c in
    ([2.0, 3.0], [1.0, -2.0, 1.5])
    taylor = vector_decay_study(RandomTaylor15(), c)
    euler = vector_decay_study(RandomEM(), c)
    @test 1.5 < taylor.strong_order < 2.3
    @test taylor.strong[end] < euler.strong[end] / 20
end

@testset "RandomTaylor15 with one noise component matches scalar noise" begin
    path = wiener_path(MersenneTwister(20261004))
    component(u, p, t, W) = -u * cos(5 * W[1])
    n = TAYLOR_STEP_COUNTS[3]
    scalar = decay_solution(RandomTaylor15(), NoiseGrid(FINE_GRID, path), n, 1.0).u
    prob = RODEProblem{false}(
        component, 1.0, (0.0, TEND), noise = NoiseGrid(FINE_GRID, [[w] for w in path])
    )
    vector = solve(prob, RandomTaylor15(), dt = TEND / n, save_everystep = true, adaptive = false).u
    @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(scalar, vector))
end

function coupled_oop(u, p, t, W)
    return [
        -u[1] + sin(u[2] + W[1]) * cos(t) + 0.3 * W[1] * W[2],
        -u[2] + cos(u[1] * W[2]) + 0.2 * t * W[1],
    ]
end
function coupled_iip(du, u, p, t, W)
    du[1] = -u[1] + sin(u[2] + W[1]) * cos(t) + 0.3 * W[1] * W[2]
    du[2] = -u[2] + cos(u[1] * W[2]) + 0.2 * t * W[1]
    return nothing
end

@testset "RandomTaylor15 in-place matches out-of-place with vector noise" begin
    rng = MersenneTwister(20261004)
    path = [[a, b] for (a, b) in zip(wiener_path(rng), wiener_path(rng))]
    u0 = [1.0, 0.5]
    for (tspan, grid, W) in (
            ((0.0, TEND), FINE_GRID, path), ((TEND, 0.0), reverse(FINE_GRID), reverse(path)),
        )
        noise = NoiseGrid(grid, W)
        dt = (tspan[2] - tspan[1]) / TAYLOR_STEP_COUNTS[3]
        outofplace = solve(
            RODEProblem{false}(coupled_oop, u0, tspan, noise = noise), RandomTaylor15(),
            dt = dt, save_everystep = true, adaptive = false
        ).u
        inplace = solve(
            RODEProblem(coupled_iip, copy(u0), tspan, noise = noise), RandomTaylor15(),
            dt = dt, save_everystep = true, adaptive = false
        ).u
        @test all(isapprox(a, b, rtol = 1.0e-9) for (a, b) in zip(outofplace, inplace))
    end
end

@testset "RandomTaylor15 vector step integrals are exact on a misaligned grid" begin
    grid = collect(range(0.0, 1.0; length = 13))
    rng = MersenneTwister(20261004)
    values = [randn(rng, 2) for _ in grid]
    path(s) = (k = min(floor(Int, 12 * s) + 1, 12); values[k] .+ (values[k + 1] .- values[k]) .* (12 * s - (k - 1)))
    integrals! = StochasticDiffEq.StochasticDiffEqRODE.path_integrals!
    for (t, dt) in ((0.0, 0.25), (0.05, 0.25), (0.03, 0.04), (1 / 12, 1 / 12), (0.5, 0.5))
        for (g, v, t0, h) in ((grid, values, t, dt), (reverse(grid), reverse(values), t + dt, -dt))
            w0 = path(t0)
            W = (t = g, W = v, dW = path(t0 + h) .- w0, curW = w0)
            I1 = zeros(2)
            I2 = zeros(2, 2)
            integrals!(I1, I2, ones(2), zeros(2), W, t0, h, w0)
            knots = sort(unique([t0; t0 + h; filter(s -> min(t0, t0 + h) < s < max(t0, t0 + h), grid)]))
            R1 = zeros(2)
            R2 = zeros(2, 2)
            for (a, b) in zip(knots[1:(end - 1)], knots[2:end])
                va, vm, vb = path(a) .- w0, path((a + b) / 2) .- w0, path(b) .- w0
                R1 .+= sign(h) * (b - a) / 6 .* (va .+ 4 .* vm .+ vb)
                R2 .+= sign(h) * (b - a) / 6 .* (va * va' .+ 4 .* (vm * vm') .+ vb * vb')
            end
            @test I1 ≈ R1
            @test I2 ≈ R2
        end
    end
end

@testset "RandomTaylor15 integrates backwards in time with vector noise" begin
    rng = MersenneTwister(20261004)
    c = [2.0, 3.0]
    paths = [wiener_path(rng), wiener_path(rng)]
    exact = cosine_decay_exact(c[1] .* paths[1] .+ c[2] .* paths[2])
    decay(u, p, t, W) = -u * cos(sum(c .* W))
    n = TAYLOR_STEP_COUNTS[3]
    noise = NoiseGrid(reverse(FINE_GRID), reverse([[a, b] for (a, b) in zip(paths...)]))
    prob = RODEProblem{false}(decay, exact[end], (TEND, 0.0), noise = noise)
    sol = solve(prob, RandomTaylor15(), dt = -TEND / n, save_everystep = true, adaptive = false)
    stride = FINE_POINTS ÷ n
    @test maximum(abs(sol.u[k + 1] - exact[end - k * stride]) for k in 0:n) < 2.0e-3
end
