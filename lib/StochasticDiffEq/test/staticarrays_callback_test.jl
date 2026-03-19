# Test for DifferentialEquations.jl#1129:
# ContinuousCallback with SDEProblem using StaticArrays (SVector)
# This should not attempt to mutate the SVector state
#
# The bug was that sde_interpolant! was being called for out-of-place problems,
# which tried to mutate immutable SVectors. The fix ensures that out-of-place
# problems use non-mutating interpolation (linear_interpolant) via the
# is_constant_cache check in change_t_via_interpolation!.

using StochasticDiffEq, StaticArrays, Test, SciMLBase

@testset "ContinuousCallback with SDEProblem{SVector} (issue #1129)" begin
    # Define out-of-place drift and diffusion functions
    function f(u::V, p, t) where {V <: StaticVector}
        V(u[2], -p)
    end
    function g(::V, p, t) where {V <: StaticVector}
        V(0.0, 0.0)
    end

    # Condition: trigger when u[1] crosses zero
    function condition(u, t, integrator)
        u[1]
    end

    u0 = SVector(1.0, -2.0)
    tspan = (0.0, 2.0)
    p = 0.0

    # Create out-of-place SDEProblem
    fn = SDEFunction{false}(f, g)
    sdeprob = SDEProblem(fn, u0, tspan, p)

    # First verify it works without callback
    sol_no_cb = solve(sdeprob, SOSRI())
    @test sol_no_cb.retcode == ReturnCode.Success

    # Now test with ContinuousCallback - this is what triggers the bug
    cb = ContinuousCallback(condition, terminate!)

    # This should not error - should use non-mutating interpolation
    # Before the fix, this would throw:
    # "setindex!(::StaticArraysCore.SVector{2, Float64}, value, ::Int) is not defined"
    sol = solve(sdeprob, SOSRI(); callback = cb)
    @test sol.retcode == ReturnCode.Terminated

    # Verify the callback actually triggered (u[1] should be near zero at the end)
    @test abs(sol.u[end][1]) < 0.1

    # Test with different SDE solvers that support callbacks
    # Use zero noise so the callback deterministically triggers
    for alg in [EM(), LambaEM(), SOSRI(), SOSRA()]
        sol_alg = solve(sdeprob, alg; callback = cb, dt = 0.001)
        @test sol_alg.retcode == ReturnCode.Terminated
    end
end

@testset "DiscreteCallback with SDEProblem{SVector}" begin
    # Also test DiscreteCallback with StaticArrays for completeness
    function f(u::V, p, t) where {V <: StaticVector}
        V(-u[1], -u[2])
    end
    function g(u::V, p, t) where {V <: StaticVector}
        V(0.1, 0.1)
    end

    u0 = SVector(1.0, 1.0)
    tspan = (0.0, 1.0)

    fn = SDEFunction{false}(f, g)
    sdeprob = SDEProblem(fn, u0, tspan)

    # Simple DiscreteCallback that terminates at t=0.5
    condition_dc(u, t, integrator) = t >= 0.5
    affect_dc!(integrator) = terminate!(integrator)
    cb_dc = DiscreteCallback(condition_dc, affect_dc!)

    sol = solve(sdeprob, SOSRI(); callback = cb_dc)
    @test sol.retcode == ReturnCode.Terminated
    @test sol.t[end] ≈ 0.5 atol = 0.01
end
