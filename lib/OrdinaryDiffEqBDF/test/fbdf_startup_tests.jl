using OrdinaryDiffEqBDF, LinearAlgebra, Test, ODEProblemLibrary
using OrdinaryDiffEqNonlinearSolve: NLFunctional
import OrdinaryDiffEqCore, LinearSolve, SciMLBase

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
scaled_decay!(du, u, p, t) = (du .= -p .* u)   # with mass matrix p·I: still y′ = -y

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

@testset "Step growth is limited by qmax" begin
    # Stiff Van der Pol starts from a tiny initial dt; without a limit the second step grew
    # by about 1e10. Growth is limited to qmax_first_step after the first step and to qmax
    # per step after that. The last step is clipped to tspan and left out.
    prob = ODEProblemLibrary.prob_ode_vanderpol_stiff
    qmax_first = OrdinaryDiffEqCore.qmax_first_step_default(FBDF())
    for qmax in (10 // 1, 5 // 1)   # 10 // 1 is the default; Rational like the qsteady defaults
        sol = solve(prob, FBDF(; qmax); reltol = 1.0e-6, abstol = 1.0e-6)
        @test sol.retcode == ReturnCode.Success
        d = diff(sol.t)[1:(end - 1)]
        growth = d[2:end] ./ d[1:(end - 1)]
        # diff(sol.t) carries the rounding of t, hence the 1e-4 margin.
        @test growth[1] <= qmax_first * (1 + 1.0e-4)
        @test maximum(growth[2:end]) <= qmax * (1 + 1.0e-4)
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

@testset "FBDF startup supports multidimensional state arrays" begin
    prob = ODEProblem(decay!, ones(2, 2, 2), (0.0, 1.0))
    for alg in (FBDF(), FBDF(nlsolve = NLFunctional()))
        sol = solve(prob, alg; reltol = 1.0e-6, abstol = 1.0e-6)
        @test sol.retcode == ReturnCode.Success
        @test size(sol.u[end]) == size(prob.u0)
        @test sol.u[end] ≈ exp(-1) .* prob.u0 rtol = 1.0e-5
    end
end

# Fail only the extra error-estimator solve, after Newton has converged.
# A failed solver is allowed to overwrite its input and output buffers.
struct FailFBDFEstimate <: LinearSolve.AbstractSolveFunction
    rhs::Base.RefValue{Any}
    failures::Base.RefValue{Int}
    persistent::Bool
end

LinearSolve.needs_concrete_A(::FailFBDFEstimate) = true

function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::FailFBDFEstimate; kwargs...)
    if cache.b === alg.rhs[] && (alg.persistent || alg.failures[] == 0)
        alg.failures[] += 1
        fill!(cache.b, NaN)
        fill!(cache.u, NaN)
        return SciMLBase.build_linear_solution(
            alg, cache.u, nothing, nothing; retcode = ReturnCode.ConvergenceFailure
        )
    end
    cache.u .= cache.A \ cache.b
    return SciMLBase.build_linear_solution(
        alg, cache.u, nothing, nothing; retcode = ReturnCode.Success
    )
end

# An FBDF integrator whose first-step estimator solve fails (once, or every time).
function init_failing_estimate(prob, persistent; kwargs...)
    linsolve = FailFBDFEstimate(Ref{Any}(nothing), Ref(0), persistent)
    integ = init(prob, FBDF(; linsolve); reltol = 1.0e-6, abstol = 1.0e-6, kwargs...)
    linsolve.rhs[] = integ.cache.terkp1_tmp   # identifies the estimator solve
    return integ, linsolve
end

@testset "FBDF rejects a failed first-step error-estimator solve" begin
    # Rescaling y' = -y must not make an inaccurate step acceptable. An unfiltered
    # residual would be smaller by 1e-9 even though the state error is unchanged.
    for scale in (1.0, 1.0e-9), persistent in (false, true)
        prob = ODEProblem(
            ODEFunction(scaled_decay!; mass_matrix = Diagonal([scale])), [1.0], (0.0, 1.0),
            scale
        )
        integ, linsolve = init_failing_estimate(prob, persistent; dt = 0.1, maxiters = 10)
        if persistent
            sol = solve!(integ)
            @test !SciMLBase.successful_retcode(sol)
            @test integ.stats.naccept == 0
            @test integ.t == prob.tspan[1]
            @test linsolve.failures[] > 1
        else
            step!(integ)
            @test linsolve.failures[] == 1
            @test 0 < integ.t < 0.1
            @test isfinite(OrdinaryDiffEqCore.get_EEst(integ))
            true_error = abs(integ.u[1] - exp(-integ.t)) / (2.0e-6)
            @test true_error < 1
        end
    end
end

@testset "FBDF falls back to the unfiltered estimate for M = I" begin
    # For M = I the residual is a state error, so a failing estimator solve must not stop the solve.
    integ, linsolve = init_failing_estimate(ODEProblem(decay!, [1.0], (0.0, 1.0)), true)
    sol = solve!(integ)
    @test sol.retcode == ReturnCode.Success
    @test linsolve.failures[] >= 1
    @test sol.u[end][1] ≈ exp(-1) rtol = 1.0e-4
end
