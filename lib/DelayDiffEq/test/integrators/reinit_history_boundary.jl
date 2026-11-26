using DelayDiffEq
using OrdinaryDiffEqTsit5
using SciMLBase: ReturnCode
using Test

# Regression for SciML/OrdinaryDiffEq.jl#4743.
#
# Exact solution u(t) = 2 + t with lag τ = 0.125. After reinit! to a later
# interval, history queries in [t0_old, t0_new) must use the user history h,
# not interpolation/extrapolation of the freshly reset solution.

const τ_4743 = 0.125
h_4743(p, t) = 2 + t
f_4743(u, h, p, t) = (h(p, t - τ_4743) + τ_4743) / (2 + t)
function f_4743!(du, u, h, p, t)
    du[1] = (h(p, t - τ_4743; idxs = 1) + τ_4743) / (2 + t)
    return nothing
end
h_4743_vec(p, t; idxs = nothing) = idxs === nothing ? [2 + t] : (2 + t)
function h_4743_vec(val, p, t; idxs = nothing)
    if idxs === nothing
        val[1] = 2 + t
    else
        val .= (2 + t)
    end
    return nothing
end

@testset "issue #4743: reinit! history boundary (oop)" begin
    prob = DDEProblem(f_4743, 1.0, h_4743, (-1.0, 0.0); constant_lags = [τ_4743])
    integrator = init(
        prob, MethodOfSteps(Tsit5()); dt = 0.2, abstol = 1.0e-8, reltol = 1.0e-8
    )
    reinit!(integrator, 2.0; t0 = 0.0, tf = 1.0)

    @test integrator.f.h(nothing, -0.1) == h_4743(nothing, -0.1) == 1.9
    @test integrator.f.h(nothing, -0.1, Val{0}) == 1.9

    sol = solve!(integrator)
    @test sol.retcode == ReturnCode.Success
    # With correct history, Tsit5 integrates u'=1 exactly up to floating-point noise.
    # Without the fix, the cached history returns 2.0 and max error is O(1e-3).
    @test maximum(abs(u - (2 + t)) for (u, t) in zip(sol.u, sol.t)) < 1.0e-12

    prob_fresh = DDEProblem(f_4743, 2.0, h_4743, (0.0, 1.0); constant_lags = [τ_4743])
    sol_fresh = solve(
        prob_fresh, MethodOfSteps(Tsit5()); dt = 0.2, abstol = 1.0e-8, reltol = 1.0e-8
    )
    @test sol_fresh.retcode == ReturnCode.Success
    @test maximum(abs(u - (2 + t)) for (u, t) in zip(sol_fresh.u, sol_fresh.t)) < 1.0e-12
end

@testset "issue #4743: reinit! history boundary (iip)" begin
    prob = DDEProblem(f_4743!, [1.0], h_4743_vec, (-1.0, 0.0); constant_lags = [τ_4743])
    integrator = init(
        prob, MethodOfSteps(Tsit5()); dt = 0.2, abstol = 1.0e-8, reltol = 1.0e-8
    )
    reinit!(integrator, [2.0]; t0 = 0.0, tf = 1.0)

    @test integrator.f.h(nothing, -0.1) == [1.9]
    val = [NaN]
    integrator.f.h(val, nothing, -0.1)
    @test only(val) == 1.9

    sol = solve!(integrator)
    @test sol.retcode == ReturnCode.Success
    @test maximum(abs(only(u) - (2 + t)) for (u, t) in zip(sol.u, sol.t)) < 1.0e-12
end

@testset "issue #4743: save_start=false / erase_sol=false" begin
    for save_start in (false, true), erase_sol in (false, true)
        prob = DDEProblem(f_4743, 1.0, h_4743, (-1.0, 0.0); constant_lags = [τ_4743])
        integrator = init(
            prob, MethodOfSteps(Tsit5());
            dt = 0.2, abstol = 1.0e-8, reltol = 1.0e-8, save_start
        )
        for _ in 1:3
            step!(integrator)
        end
        reinit!(integrator, 2.0; t0 = 0.0, tf = 1.0, erase_sol)

        @test integrator.f.h(nothing, -0.1) == 1.9

        sol = solve!(integrator)
        @test sol.retcode == ReturnCode.Success
        # With erase_sol=false, earlier interval points remain; check the new span.
        idx = findall(t -> t >= 0, sol.t)
        @test !isempty(idx)
        @test maximum(abs(sol.u[i] - (2 + sol.t[i])) for i in idx) < 1.0e-12
    end
end

@testset "issue #4743: backwards reinit! history boundary" begin
    # Exact solution u(t)=2+t integrated backwards. Lag looks forward in clock
    # time (negative constant_lags, matching test/interface/backwards.jl).
    f_back(u, h, p, t) = (h(p, t + τ_4743) - τ_4743) / (2 + t)
    prob = DDEProblem(f_back, 3.0, h_4743, (1.0, 0.0); constant_lags = [-τ_4743])
    integrator = init(
        prob, MethodOfSteps(Tsit5()); dt = 0.2, abstol = 1.0e-8, reltol = 1.0e-8
    )
    reinit!(integrator, 2.0; t0 = 0.0, tf = -1.0)

    # t = 0.1 is past the new start in clock time, i.e. before the start along tdir.
    @test integrator.tdir < 0
    @test integrator.f.h(nothing, 0.1) == h_4743(nothing, 0.1) == 2.1

    sol = solve!(integrator)
    @test sol.retcode == ReturnCode.Success
    @test maximum(abs(u - (2 + t)) for (u, t) in zip(sol.u, sol.t)) < 1.0e-12
end
