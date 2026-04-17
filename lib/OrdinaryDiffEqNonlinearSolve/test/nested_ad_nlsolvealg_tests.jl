using Test
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using DiffEqBase: DAEProblem, ODEProblem, ReturnCode
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using ForwardDiff

# Regression test for "No matching function wrapper was found!" errors raised
# during `init` when an implicit ODE solver using `NonlinearSolveAlg` as its
# inner nonlinear solver is invoked inside an outer ForwardDiff.Dual layer
# (i.e. the user is differentiating the full ODE solve).
#
# Before the fix, `build_nlsolver` constructed the inner NonlinearProblem with
# an unspecialized `NonlinearFunction`, which SciMLBase wrapped in a
# FunctionWrappersWrapper covering only a fixed set of Dual-tag branches
# (Float64 + NonlinearSolveBase.NonlinearSolveTag). When the ODE solver was
# itself called with a `Dual{DiffEqBase.OrdinaryDiffEqTag}` input, no matching
# wrapper branch existed and FunctionWrappersWrappers threw at dispatch time.

@testset "Nested ForwardDiff over implicit solver with NonlinearSolveAlg" begin
    function lin!(du, u, p, t)
        @. du = -p * u
        return nothing
    end

    u0 = [1.0, 2.0, 3.0]
    tspan = (0.0, 1.0)

    # Gradient of the final-time sum of state with respect to a 3-element
    # parameter vector. Each of the three parameter components produces a
    # distinct ForwardDiff chunk direction, so the chunk size here is 3 and
    # u/du in `init` become `Vector{Dual{OrdinaryDiffEqTag, Float64, 3}}`.
    function make_loss(alg)
        return function (p)
            prob = ODEProblem(lin!, u0, tspan, p)
            sol = solve(
                prob, alg; save_everystep = false, reltol = 1.0e-8,
                abstol = 1.0e-8
            )
            return sum(sol.u[end])
        end
    end

    for alg in (
            ImplicitEuler(nlsolve = NonlinearSolveAlg()),
            TRBDF2(nlsolve = NonlinearSolveAlg()),
        )
        loss = make_loss(alg)
        # The @test_nowarn just guards against a regression — the important
        # thing is that this does not throw
        # `ArgumentError: No matching function wrapper was found!`.
        g = ForwardDiff.gradient(loss, [1.0, 1.0, 1.0])
        @test length(g) == 3
        @test all(isfinite, g)
    end
end

@testset "DFBDF DAE init with internal ForwardDiff Jacobian" begin
    # Companion regression test: the DAE initialization path in
    # `initialize_dae.jl` also builds NonlinearProblems whose inner f got
    # FunctionWrappersWrapped on only Float64 / NonlinearSolveTag branches,
    # which fails when OrdinaryDiffEq's own ForwardDiff chunk Jacobian
    # (tag = `DiffEqBase.OrdinaryDiffEqTag`) calls back into that wrapper.
    function lorenz_dae!(out, du, u, p, t)
        out[1] = 10.0 * (u[2] - u[1]) - du[1]
        out[2] = u[1] * (28.0 - u[3]) - u[2] - du[2]
        out[3] = u[1] * u[2] - (8 / 3) * u[3] - du[3]
        return nothing
    end

    u0 = [1.0, 0.0, 0.0]
    du0 = [0.0, 0.0, 0.0]
    prob = DAEProblem(
        lorenz_dae!, du0, u0, (0.0, 1.0);
        differential_vars = [true, true, true]
    )
    sol = solve(prob, DFBDF(); reltol = 1.0e-6, abstol = 1.0e-6,
        initializealg = BrownFullBasicInit())
    @test sol.retcode == ReturnCode.Success
end
