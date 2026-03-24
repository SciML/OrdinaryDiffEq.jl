using OrdinaryDiffEq, ForwardDiff, FiniteDiff, SciMLSensitivity, Zygote, Test

# Regression test for ForwardDiff AD through ContinuousCallback root-finding.
# This catches issues like https://github.com/SciML/DiffEqBase.jl/issues/1275
# where the bracketing solver (e.g. ModAB) fails to handle ForwardDiff.Dual numbers
# correctly during event detection, producing NaN derivatives.

@testset "ForwardDiff AD through ContinuousCallback" begin
    function f(du, u, p, t)
        du[1] = u[2]
        du[2] = -p[1]
    end

    # Callback triggers when u[1] crosses zero (bouncing ball condition)
    # and modifies the state (reverses velocity with coefficient of restitution)
    function condition(u, t, integrator)
        u[1]
    end

    function affect!(integrator)
        integrator.u[2] = -integrator.p[2] * integrator.u[2]
    end

    cb = ContinuousCallback(condition, affect!)

    function solve_and_extract(p)
        u0 = [1.0, 0.0]
        tspan = (0.0, 3.0)
        prob = ODEProblem(f, u0, tspan, p)
        sol = solve(prob, Tsit5(), callback = cb, abstol = 1.0e-12, reltol = 1.0e-12)
        return [sol.u[end][1], sol.u[end][2]]
    end

    # Test over a range of parameter values to cover different callback timings
    @testset "p1 = $p1" for p1 in 0.5:0.5:5.0
        p = [p1, 0.8]
        dijac = ForwardDiff.jacobian(solve_and_extract, p)
        findiff = FiniteDiff.finite_difference_jacobian(solve_and_extract, p)
        @test all(isfinite, dijac)
        @test dijac ≈ findiff rtol = 1.0e-5
    end
end

# Test with mutable closure (from OrdinaryDiffEq AD tests pattern)
@testset "ForwardDiff AD through ContinuousCallback with mutable closure" begin
    function f2(du, u, p, t)
        du[1] = -p[1] * u[1]
        du[2] = p[1] * u[1] - u[2]
    end

    called = Ref(false)

    function condition2(u, t, integrator)
        u[1] - 0.5
    end

    function affect2!(integrator)
        called[] = true
        integrator.p[2] = zero(integrator.p[2])
    end

    cb2 = ContinuousCallback(condition2, affect2!)

    function solve_with_closure(x)
        called[] = false
        u0 = [1.0, 0.0]
        tspan = (0.0, 5.0)
        prob = ODEProblem(f2, u0, tspan, copy(x))
        sol = solve(prob, Tsit5(), callback = cb2, abstol = 1.0e-12, reltol = 1.0e-12)
        return [sol.u[end][1], sol.u[end][2]]
    end

    @testset "x = $x" for x in 0.1:0.1:1.0
        p = [x, 1.0]
        dijac = ForwardDiff.jacobian(solve_with_closure, p)
        findiff = FiniteDiff.finite_difference_jacobian(solve_with_closure, p)
        @test all(isfinite, dijac)
        @test dijac ≈ findiff rtol = 1.0e-5
    end
end

# Regression test for Zygote/ReverseDiff AD through ContinuousCallback.
# ReverseDiff.TrackedReal cannot be converted to Float64, so storing the
# callback residual in integrator.last_event_error requires unwrapping
# via `value()`.
@testset "Zygote AD through ContinuousCallback" begin
    function f_zyg(du, u, p, t)
        du[1] = u[2]
        du[2] = -p[1]
    end

    function condition_zyg(u, t, integrator)
        u[1]
    end

    function affect_zyg!(integrator)
        integrator.u[2] = -integrator.p[2] * integrator.u[2]
    end

    cb_zyg = ContinuousCallback(condition_zyg, affect_zyg!)

    function loss_zyg(p)
        prob = ODEProblem(f_zyg, [1.0, 0.0], (0.0, 1.0), p)
        sol = solve(prob, Tsit5(), callback = cb_zyg,
            abstol = 1e-14, reltol = 1e-14, save_everystep = false)
        return sum(sol.u[end])
    end

    p = [9.8, 0.8]
    grad = Zygote.gradient(loss_zyg, p)[1]
    findiff_grad = FiniteDiff.finite_difference_gradient(loss_zyg, p)
    @test all(isfinite, grad)
    @test grad ≈ findiff_grad rtol = 1e-3
end
