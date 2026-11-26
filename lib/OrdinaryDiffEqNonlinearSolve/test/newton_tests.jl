using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    @test sol2.stats.nf <= sol1.stats.nf + 20
end

@testset "in-place line search preserves a resolved ascent and local trial slopes" begin
    δ = eps() / 4
    function rhs!(du, u, p, t)
        a = -u[2] / δ
        a <= 1.125 || throw(DomainError(a))
        du[1] = 1
        du[2] = u[2] + (1 + a - 1.5a^2)
        return nothing
    end
    function jac!(J, u, p, t)
        fill!(J, 0)
        J[2, 2] = 1 + inv(δ)
        return nothing
    end
    prob = ODEProblem(ODEFunction(rhs!; jac = jac!), zeros(2), (0.0, 1.0))
    function predictor!(upred, uprev, p, t, dt)
        upred[1] = 1
        upred[2] = 0
        return nothing
    end
    function solve_with(linesearch)
        return solve(
            prob, ImplicitEuler(nlsolve = NLNewton(predictor = predictor!, relax = linesearch));
            dt = 1.0, adaptive = false, abstol = 0.1, reltol = 0.0,
            internalnorm = (x, t) -> maximum(abs, x)
        )
    end

    @testset "derivative callbacks" begin
        function check_slopes(ϕ, dϕ, ϕdϕ, a, ϕ0, dϕ0)
            @test dϕ0 ≈ 10 atol = 32sqrt(eps())
            @test dϕ(0.0) == dϕ0
            @test ϕdϕ(0.0) == (ϕ0, dϕ0)
            value, slope = ϕdϕ(a)
            @test value ≈ 5
            @test slope ≈ -20 atol = 64sqrt(eps())
            @test dϕ(a) ≈ -20 atol = 64sqrt(eps())
            return a, value
        end
        solve_with(check_slopes)
    end
    if Base.pkgversion(LineSearches) >= v"7.8.2"
        @test_throws LineSearches.LineSearchException solve_with(BackTracking())
    end
end

@testset "in-place line search recovers a roundoff-limited slope" begin
    δ = 2.0^-28
    function rhs!(du, u, p, t)
        a = (1 - u[1]) / δ
        du[1] = u[1] - a
        return nothing
    end
    jac!(J, u, p, t) = (J[1, 1] = 1 + inv(δ); nothing)
    prob = ODEProblem(ODEFunction(rhs!; jac = jac!), [1.0], (0.0, 1.0))
    function check_slopes(ϕ, dϕ, ϕdϕ, a, ϕ0, dϕ0)
        @test ϕ(sqrt(eps())) == ϕ0
        @test dϕ0 ≈ -10
        @test dϕ(0.0) == dϕ0
        @test ϕdϕ(0.0) == (ϕ0, dϕ0)
        value = ϕ(a)
        slope = (ϕ(a + sqrt(eps())) - value) / sqrt(eps())
        @test dϕ(a) == slope
        @test ϕdϕ(a) == (value, slope)
        return BackTracking()(ϕ, dϕ, ϕdϕ, a, ϕ0, dϕ0)
    end
    sol = solve(
        prob, ImplicitEuler(nlsolve = NLNewton(relax = check_slopes));
        dt = 1.0, adaptive = false, abstol = 0.1, reltol = 0.0
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end] == [1 - δ]
end

@testset "fixed-step ImplicitEuler respects solve-level adaptivity (#4905)" begin
    rhs(u, p, t) = -8u
    prob = ODEProblem(rhs, 1.0, (0.0, 1.5))
    dt = 1 / 8
    integ = init(
        prob, ImplicitEuler(); adaptive = false, dt = dt,
        abstol = 1.0e-3, reltol = 1.0e-3
    )

    # Backward Euler gives u[n+1] = u[n] / (1 + 8dt) = u[n] / 2.
    expected = 1.0
    for n in 1:12
        previous_t = integ.t
        previous_iterations = integ.stats.nnonliniter
        step!(integ)

        # Stop after the first rejected step so a failed solve cannot leave the
        # test repeating at the same time.
        if !integ.accept_step
            @test integ.accept_step
            break
        end

        @test integ.t > previous_t
        @test integ.t ≈ n * dt

        iterations = integ.stats.nnonliniter - previous_iterations
        @test iterations >= 2

        expected /= 2
        @test integ.u ≈ expected atol = 1.0e-12
    end
    @test integ.stats.nreject == 0
end
