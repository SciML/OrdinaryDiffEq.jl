using OrdinaryDiffEqBDF, LinearAlgebra, Test
using OrdinaryDiffEqNonlinearSolve: NLFunctional
import OrdinaryDiffEqCore

# Start-up of FBDF (https://github.com/SciML/OrdinaryDiffEq.jl/issues/5003): the first step
# after a (re)start has no history, and order-1 steps used to be sized from ‖h y′‖.

# Capacitor discharging into an RL load, the current I as an algebraic unknown:
# V′ = -I/C, λ′ = V - R I, 0 = λ - L I, starting from rest (I = λ = 0).
const C_, L_, R_, V0_ = 2.0e-3, 1.0e-5, 0.05, 2000.0
const RLC_A = [0 0 -1 / C_; 1 0 -R_; 0 1 -L_]
const RLC_M = Diagonal([1.0, 1.0, 0.0])
const RLC_U0 = [V0_, 0.0, 0.0]
rlc!(du, u, p, t) = mul!(du, RLC_A, u)
# Exact solution: (V, λ) solve a linear ODE, I = λ / L.
function rlc_exact(t)
    x = exp([0 -1 / (C_ * L_); 1 -R_ / L_] * t) * [V0_, 0.0]
    return [x; x[2] / L_]
end
rlc_prob(t0) = ODEProblem(ODEFunction(rlc!; mass_matrix = RLC_M), RLC_U0, (t0, t0 + 1.0e-4))
rlc_solve(prob) = solve(prob, FBDF(); reltol = 1.0e-6, abstol = [1.0e-6, 1.0e-12, 1.0e-9])

decay(u, p, t) = -u
decay!(du, u, p, t) = (du .= -u)

# Weighted RMS norm of the true local error of one backward-Euler step from u0, with
# abstol = reltol = 1: an independent reference for EEst.
function true_EEst(M, A, u0, exact, h)
    u = (M - h * A) \ (M * u0)
    return sqrt(sum(abs2, (u - exact(h)) ./ (1 .+ max.(abs.(u0), abs.(u)))) / length(u))
end

@testset "First-step error estimate is the backward-Euler local error" begin
    for h in (1.0e-7, 1.0e-9)
        # Huge tolerances so the first attempt is accepted at the given dt.
        integ = init(rlc_prob(0.0), FBDF(); dt = h, reltol = 1.0, abstol = 1.0)
        step!(integ)
        @test integ.t == h
        @test OrdinaryDiffEqCore.get_EEst(integ) ≈
            true_EEst(RLC_M, RLC_A, RLC_U0, rlc_exact, h) rtol = 0.02
    end
    # y′ = -y: in place with Newton (filtered), in place with functional iteration and
    # out of place (both unfiltered).
    expected = true_EEst(I, -ones(1, 1), [1.0], t -> [exp(-t)], 1.0e-3)
    for (f, alg) in ((decay!, FBDF()), (decay!, FBDF(nlsolve = NLFunctional())), (decay, FBDF()))
        integ = init(ODEProblem(f, [1.0], (0.0, 1.0)), alg; dt = 1.0e-3, reltol = 1.0, abstol = 1.0)
        step!(integ)
        @test OrdinaryDiffEqCore.get_EEst(integ) ≈ expected rtol = 0.02
    end
end

sol0 = rlc_solve(rlc_prob(0.0))
dts0 = diff(sol0.t)

@testset "Start-up does not depend on the absolute start time" begin
    for t0 in (1.0e-8, 1.0e-3, 1.0)
        sol = rlc_solve(rlc_prob(t0))
        # Only rounding of t differs (eps(1.0) / dt ≈ 1e-6 at t0 = 1), so the start-up
        # must be the same to well below the factor-of-1e3 to 1e8 differences seen before.
        @test diff(sol.t)[1:10] ≈ dts0[1:10] rtol = 1.0e-3
        @test abs(sol.stats.naccept - sol0.stats.naccept) <= 3
        @test sol.u[end] ≈ rlc_exact(1.0e-4) rtol = 1.0e-5
    end
end

@testset "No start-up collapse on the mass-matrix DAE" begin
    # The first step is limited by backward Euler's local error (h ≈ 3e-10 s here); the
    # steps after it may not collapse below it or jump by orders of magnitude.
    @test minimum(dts0[1:(end - 1)]) >= dts0[1] / 2
    @test maximum(dts0[2:end] ./ dts0[1:(end - 1)]) < 20
    @test sol0.stats.naccept < 60
    @test sol0.u[end] ≈ rlc_exact(1.0e-4) rtol = 1.0e-5
    # Out of place, a mass-matrix problem keeps the history-based first-step estimate.
    oop = ODEProblem(ODEFunction((u, p, t) -> RLC_A * u; mass_matrix = RLC_M), RLC_U0, (0.0, 1.0e-4))
    sol = solve(oop, FBDF(); reltol = 1.0e-6, abstol = [1.0e-6, 1.0e-12, 1.0e-9])
    @test sol.retcode == ReturnCode.Success
    @test sol.u[end] ≈ rlc_exact(1.0e-4) rtol = 1.0e-4
end

@testset "Order-1 steps are sized from the local error" begin
    # Restricted to order 1, the step count must scale like tol^(-1/2) (local error h²),
    # not like 1/tol, which sizing from ‖h y′‖ gives (1.3e6 steps at tol = 1e-6).
    prob = ODEProblem(decay, [1.0], (0.0, 1.0))
    sol(tol) = solve(prob, FBDF(max_order = Val(1)); reltol = tol, abstol = tol, maxiters = 10_000)
    s4, s6 = sol(1.0e-4), sol(1.0e-6)
    @test s6.retcode == ReturnCode.Success
    @test 5 < s6.stats.naccept / s4.stats.naccept < 20
end
