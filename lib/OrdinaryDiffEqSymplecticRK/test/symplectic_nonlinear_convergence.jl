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
