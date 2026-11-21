# Nonlinear, time-dependent problems: the harmonic oscillator in
# symplectic_convergence.jl cannot detect wrong stage times or a drift/kick swap.
using OrdinaryDiffEqSymplecticRK, OrdinaryDiffEqTsit5, Test

const ALGS = (
    (PseudoVerletLeapfrog(), 2), (McAte2(), 2), (Ruth3(), 3), (McAte3(), 3),
    (CandyRoz4(), 4), (McAte4(), 4), (CalvoSanz4(), 4), (McAte42(), 4), (McAte5(), 5),
)

# Runge-Kutta-Nystrom methods: their order needs u' = v (McLachlan & Atela 1992, 3.2).
const RKN_ONLY = (McAte4, CalvoSanz4, McAte5)

# (name, kick f1, drift f2, drift depends on t)
const PROBLEMS = (
    ("pendulum", (v, u, t) -> -sin.(u), (v, u, t) -> v, false),
    ("forced pendulum", (v, u, t) -> -sin.(u) .+ cos(2t + 0.7) / 2, (v, u, t) -> v, false),
    ("time-dependent drift", (v, u, t) -> -sin.(u), (v, u, t) -> v .+ sin(3t) / 2, true),
)

function make_problem(f1, f2, iip, v0, u0, tspan)
    if iip
        return DynamicalODEProblem(
            (dv, v, u, p, t) -> (dv .= f1(v, u, t)),
            (du, v, u, p, t) -> (du .= f2(v, u, t)), v0, u0, tspan
        )
    else
        return DynamicalODEProblem(
            (v, u, p, t) -> f1(v, u, t), (v, u, p, t) -> f2(v, u, t), v0, u0, tspan
        )
    end
end

@testset "$name, iip = $iip, $(nameof(typeof(alg)))" for (name, f1, f2, drift) in PROBLEMS,
        iip in (false, true), (alg, order) in ALGS

    drift && typeof(alg) in RKN_ONLY && continue
    prob = make_problem(f1, f2, iip, [0.3], [1.1], (0.0, 1.0))
    ref = solve(prob, Tsit5(); abstol = 1.0e-14, reltol = 1.0e-14).u[end]
    dts = order <= 3 ? (1 / 32, 1 / 64) : (1 / 8, 1 / 16)
    e = [
        maximum(abs.(solve(prob, alg; dt, adaptive = false).u[end] .- ref))
            for dt in dts
    ]
    @test log2(e[1] / e[2]) ≈ order atol = 0.3
end

@testset "endpoint time, iip = $iip, $(nameof(typeof(alg)))" for iip in (false, true),
        (alg, _) in ALGS

    seen = Float64[]
    f1(v, u, t) = (push!(seen, t); zero(v))
    f2(v, u, t) = v
    prob = make_problem(f1, f2, iip, [0.0], [0.0], (1.0, 1.2))
    solve(prob, alg; dt = 0.2, adaptive = false)
    @test seen[end] == 1.2
end

# Stage-clock invariant: kick i at t + (b[1] + ... + b[i]) * dt, the last at t + dt.
@testset "kick times, iip = $iip, $(nameof(typeof(alg)))" for iip in (false, true),
        (alg, _) in ALGS

    tab = OrdinaryDiffEqSymplecticRK._symplectic_tableau(alg, Float64, Float64)
    N = length(tab.a)
    expected = [cumsum(collect(tab.b))[1:(N - 1)]; 1.0]
    iszero(tab.a[N]) || push!(expected, 1.0)
    seen = Float64[]
    f1(v, u, t) = (push!(seen, t); zero(v))
    f2(v, u, t) = v
    integrator = init(
        make_problem(f1, f2, iip, [0.0], [0.0], (0.0, 1.0)), alg; dt = 1.0, adaptive = false
    )
    empty!(seen)
    step!(integrator)
    @test seen ≈ expected atol = 1.0e-14
    @test seen[end] == 1.0
end

# Ruth3 and CandyRoz4 evaluate the force outside the step (documented); the
# suggested alternatives stay inside it.
@testset "callback times at the boundaries, iip = $iip" for iip in (false, true)
    function callback_times(alg)
        kicks, drifts = Float64[], Float64[]
        f1(v, u, t) = (push!(kicks, t); -sin.(u))
        f2(v, u, t) = (push!(drifts, t); v)
        solve(
            make_problem(f1, f2, iip, [0.3], [1.1], (0.0, 1.0)), alg;
            dt = 0.25, adaptive = false
        )
        return kicks, drifts
    end
    kicks, _ = callback_times(Ruth3())
    @test maximum(kicks) ≈ 0.75 + 0.25 * 25 / 24
    @test minimum(kicks) >= 0
    kicks, _ = callback_times(CandyRoz4())
    @test minimum(kicks) ≈ -0.25 * 0.35120719195965777
    @test maximum(kicks) ≈ 0.75 + 0.25 * 1.3512071919596578
    for alg in (PseudoVerletLeapfrog(), McAte2(), McAte3(), KahanLi6(), KahanLi8(), SofSpa10())
        kicks, drifts = callback_times(alg)
        @test all(t -> 0 <= t <= 1, kicks)
        @test all(t -> 0 <= t <= 1, drifts)
    end
end

# Stage times are formed from t, not accumulated, so rounding cannot push an
# in-range stage outside a tiny step.
@testset "tiny step stays inside tspan, iip = $iip" for iip in (false, true)
    tspan = (1.0, nextfloat(1.0, 8))
    inside(t) = tspan[1] <= t <= tspan[2] || throw(DomainError(t))
    f1(v, u, t) = (inside(t); zero(v))
    f2(v, u, t) = (inside(t); v)
    for alg in (PseudoVerletLeapfrog(), McAte2(), McAte3(), KahanLi6(), KahanLi8(), SofSpa10())
        sol = solve(
            make_problem(f1, f2, iip, [0.0], [1.0], tspan), alg;
            dt = 8 * eps(1.0), adaptive = false
        )
        @test sol.t[end] == tspan[2]
    end
end

# No evaluation past tspan when t + dt rounds beyond it or the coefficient sum rounds above 1.
const INSIDE_ALGS = (
    SymplecticEuler(), VelocityVerlet(), VerletLeapfrog(), LeapfrogDriftKickDrift(),
    PseudoVerletLeapfrog(), McAte2(), McAte3(), KahanLi6(), KahanLi8(), SofSpa10(),
)

function guarded_solve(alg, iip, tspan, dt, T)
    lo, hi = minmax(tspan...)
    inside(t) = lo <= t <= hi || throw(DomainError(t))
    f1(v, u, t) = (inside(t); -sin.(u))
    f2(v, u, t) = (inside(t); v)
    return solve(
        make_problem(f1, f2, iip, T[0.3], T[1.1], tspan), alg; dt, adaptive = false
    )
end

@testset "non-binary endpoint, iip = $iip, $(nameof(typeof(alg)))" for iip in (false, true),
        alg in INSIDE_ALGS

    sol = guarded_solve(alg, iip, (-1.0, 0.1), 1.1, Float64)
    @test sol.t[end] == 0.1
    sol = guarded_solve(alg, iip, (0.1, -1.0), 1.1, Float64)
    @test sol.t[end] == -1.0
end

@testset "Float32 endpoint, iip = $iip, $(nameof(typeof(alg)))" for iip in (false, true),
        alg in INSIDE_ALGS

    sol = guarded_solve(alg, iip, (0.0f0, 1.0f0), 1.0f0, Float32)
    @test sol.t[end] == 1.0f0
    sol = guarded_solve(alg, iip, (0.0f0, 1.0f0), 0.1f0, Float32)
    @test sol.t[end] == 1.0f0
end
