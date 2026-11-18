# Order on nonlinear problems, with and without a time-dependent force.
# The harmonic oscillator in symplectic_convergence.jl is linear and autonomous, so it
# does not detect a wrong time for the kicks (#5014) or a drift/kick swap that only
# matters when the force is nonlinear (McAte4, CalvoSanz4, McAte5, McAte42 #85).
using OrdinaryDiffEqSymplecticRK, OrdinaryDiffEqTsit5, Test

pendulum(v, u, p, t) = -sin.(u)
forced_pendulum(v, u, p, t) = -sin.(u) .+ cos(2t + 0.7) / 2

const ALGS = (
    (PseudoVerletLeapfrog(), 2), (McAte2(), 2), (Ruth3(), 3), (McAte3(), 3),
    (CandyRoz4(), 4), (McAte4(), 4), (CalvoSanz4(), 4), (McAte42(), 4), (McAte5(), 5),
)

@testset "$name, $(nameof(typeof(alg)))" for (name, force) in
        (("pendulum", pendulum), ("forced pendulum", forced_pendulum)), (alg, order) in ALGS
    prob = SecondOrderODEProblem(force, [0.3], [1.1], (0.0, 1.0))
    ref = solve(prob, Tsit5(); abstol = 1.0e-14, reltol = 1.0e-14).u[end]
    dts = order <= 3 ? (1 / 32, 1 / 64) : (1 / 8, 1 / 16)
    e = [maximum(abs.(solve(prob, alg; dt, adaptive = false).u[end] .- ref)) for dt in dts]
    @test log2(e[1] / e[2]) ≈ order atol = 0.3
end

# Time-dependent drift: u' = v + sin(3t)/2, v' = -sin(u). The drift is evaluated at
# the velocity's time and the kick at the position's time. McAte4, CalvoSanz4 and
# McAte5 are Runge-Kutta-Nystrom methods (McLachlan & Atela 1992, section 3.2): their
# order needs u' = v, so they are left out here.
const RKN_ONLY = (McAte4, CalvoSanz4, McAte5)
const DRIFT_ALGS = filter(x -> !(typeof(x[1]) in RKN_ONLY), ALGS)
drift_f1(v, u, p, t) = -sin.(u)
drift_f2(v, u, p, t) = v .+ sin(3t) / 2

@testset "time-dependent drift, $(nameof(typeof(alg)))" for (alg, order) in DRIFT_ALGS
    prob = DynamicalODEProblem(drift_f1, drift_f2, [0.3], [1.1], (0.0, 1.0))
    ref = solve(prob, Tsit5(); abstol = 1.0e-14, reltol = 1.0e-14).u[end]
    dts = order <= 3 ? (1 / 32, 1 / 64) : (1 / 8, 1 / 16)
    e = [
        maximum(abs.(solve(prob, alg; dt, adaptive = false).u[end] .- ref))
            for dt in dts
    ]
    @test log2(e[1] / e[2]) ≈ order atol = 0.3
end

# The step ends exactly at t + dt, so the cached derivatives are evaluated there.
@testset "endpoint time, $(nameof(typeof(alg)))" for (alg, _) in ALGS
    seen = Float64[]
    f1(v, u, p, t) = (push!(seen, t); zero(v))
    f2(v, u, p, t) = v
    prob = DynamicalODEProblem(f1, f2, [0.0], [0.0], (1.0, 1.2))
    solve(prob, alg; dt = 0.2, adaptive = false)
    @test seen[end] == 1.2
end
