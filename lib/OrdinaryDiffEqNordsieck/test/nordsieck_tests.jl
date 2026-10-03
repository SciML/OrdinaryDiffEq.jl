using OrdinaryDiffEqNordsieck, DiffEqDevTools, Test, LinearAlgebra
import ODEProblemLibrary: prob_ode_bigfloatlinear,
    prob_ode_bigfloat2Dlinear,
    prob_ode_linear, prob_ode_2Dlinear

probArr = [prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear]
testTol = 0.25
dts = 1 .// (2 .^ (10:-1:5))

@testset "Nordsieck Convergence Tests" begin
    for i in eachindex(probArr)
        sim = test_convergence(dts, probArr[i], AN5())
        @test sim.𝒪est[:final] ≈ 5 atol = testTol
        @test sim.𝒪est[:l2] ≈ 5 atol = testTol
        @test sim.𝒪est[:l∞] ≈ 5 atol = testTol
    end
end

probArr = [
    prob_ode_linear,
    prob_ode_2Dlinear,
]
@testset "Nordsieck Adaptivity Tests: AN5" begin
    for i in eachindex(probArr)
        prob = probArr[i]
        sol = solve(prob, AN5(), reltol = 1.0e-6)
        @test length(sol.t) < 11
        @test SciMLBase.successful_retcode(sol)
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test exact ≈ sol.u[end] atol = 1.0e-5
    end
end

@testset "Nordsieck Adaptivity Tests: JVODE" begin
    for i in eachindex(probArr),
            sol in [JVODE_Adams(), JVODE_BDF()]

        prob = probArr[i]
        sol = solve(prob, sol, reltol = 1.0e-4, abstol = 1.0e-7)
        @test length(sol.t) < 22
        @test SciMLBase.successful_retcode(sol)
        exact = prob.f.analytic(prob.u0, prob.p, prob.tspan[end])
        @test norm(exact - sol.u[end], Inf) < 3.0e-3
    end
end

# Coupled linear resize: AN5Cache had no `full_cache`, so `resize!(integrator, n)`
# threw MethodError. History after a resize is invalid, so the method restarts at
# order 1. JVODE is omitted: its `post_newton_controller!` 2-arg/3-arg recursion
# (StackOverflowError) is a separate bug, not this cache-resize path.
const A_resize = [-1.0 0.35 0.1; 0.2 -1.4 0.25; 0.15 0.3 -1.2]
const U0_resize = [1.0, -0.5]
const UNEW_resize = 0.3
const TG_resize = 0.5
const TS_resize = 0.65
const TF_resize = 0.8

f_resize!(du, u, p, t) = (mul!(du, view(A_resize, 1:length(u), 1:length(u)), u); nothing)

function exact_resize(t)
    u = exp(A_resize[1:2, 1:2] * min(t, TG_resize)) * U0_resize
    t <= TG_resize && return u
    u3 = exp(A_resize * (min(t, TS_resize) - TG_resize)) * vcat(u, UNEW_resize)
    t <= TS_resize && return u3
    return exp(A_resize[1:2, 1:2] * (t - TS_resize)) * u3[1:2]
end

@testset "AN5 state resizing" begin
    grew = Ref(false)
    shrank = Ref(false)
    restart_orders = Int[]
    function affect!(integrator)
        n = grew[] ? 2 : 3
        resize!(integrator, n)
        n == 3 && (integrator.u[3] = UNEW_resize)
        # New `uprev` slots are uninitialized (#4722); a Nordsieck restart
        # reads `uprev` as the step start, so it must match the resized `u`.
        copyto!(integrator.uprev, integrator.u)
        push!(restart_orders, integrator.cache.order)
        grew[] ? (shrank[] = true) : (grew[] = true)
        return nothing
    end
    cb = DiscreteCallback(
        (u, t, integrator) -> (!grew[] && t >= TG_resize) || (!shrank[] && t >= TS_resize),
        affect!;
        save_positions = (false, false)
    )
    sol = solve(
        ODEProblem(f_resize!, copy(U0_resize), (0.0, TF_resize)), AN5();
        callback = cb, tstops = [TG_resize, TS_resize], reltol = 1.0e-10, abstol = 1.0e-12
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol(0.58) ≈ exact_resize(0.58) rtol = 1.0e-7 atol = 1.0e-10
    @test sol.u[end] ≈ exact_resize(TF_resize) rtol = 1.0e-7 atol = 1.0e-10
    @test restart_orders == [1, 1]
end
