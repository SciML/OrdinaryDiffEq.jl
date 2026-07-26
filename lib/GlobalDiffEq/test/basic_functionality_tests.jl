using GlobalDiffEq, OrdinaryDiffEqTsit5, LinearAlgebra
using OrdinaryDiffEqSSPRK
using Test

@testset "Basic functionality" begin
    l = 1.0                             # length [m]
    m = 1.0                             # mass[m]
    g = 9.81                            # gravitational acceleration [m/s²]

    function pendulum!(du, u, p, t)
        du[1] = u[2]                    # θ'(t) = ω(t)
        return du[2] = -3g / (2l) * sin(u[1]) + 3 / (m * l^2) * p(t) # ω'(t) = -3g/(2l) sin θ(t) + 3/(ml^2)M(t)
    end

    θ₀ = 0.01                           # initial angular deflection [rad]
    ω₀ = 0.0                            # initial angular velocity [rad/s]
    u₀ = [θ₀, ω₀]                       # initial state vector
    tspan = (0.0, 10.0)                  # time interval

    M = t -> 0.1sin(t)                    # external torque [Nm]

    prob = ODEProblem(pendulum!, u₀, tspan, M)

    v0 = solve(prob, Tsit5(), dt = 0.1, reltol = 1.0e-12, abstol = 0).(1:10)
    ve = solve(
        prob, GlobalRichardson(SSPRK33()), dt = 0.2, reltol = 1.0e-12, abstol = 0
    ).(1:10)
    @test norm(ve - v0) / norm(v0) < 1.0e-10
end
