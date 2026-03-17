using DelayDiffEq
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner
using OrdinaryDiffEqBDF
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using LinearAlgebra
using Test

@testset "SIR" begin
    function sir_dde!(du, u, h, p, t)
        S, I, _ = u
        γ, τ = p
        infection = γ * I * S
        Sd, Id, _ = h(p, t - τ)
        recovery = γ * Id * Sd
        @inbounds begin
            du[1] = -infection
            du[2] = infection - recovery
            du[3] = recovery
        end
        nothing
    end
    u0 = [0.99, 0.01, 0.0]
    sir_history(p, t) = [1.0, 0.0, 0.0]
    tspan = (0.0, 40.0)
    p = (γ = 0.5, τ = 4.0)

    prob_dde = DDEProblem(sir_dde!, u0, sir_history, tspan, p; constant_lags = (p.τ,))
    sol_dde = TestSolution(
        solve(
            prob_dde, MethodOfSteps(Vern9()); reltol = 1.0e-14,
            abstol = 1.0e-14
        )
    )

    function sir_ddae!(du, u, h, p, t)
        S, I, R = u
        γ, τ = p
        infection = γ * I * S
        Sd, Id, _ = h(p, t - τ)
        recovery = γ * Id * Sd
        @inbounds begin
            du[1] = -infection
            du[2] = infection - recovery
            du[3] = S + I + R - 1
        end
        nothing
    end

    prob_ddae = DDEProblem(
        DDEFunction{true}(
            sir_ddae!;
            mass_matrix = Diagonal([1.0, 1.0, 0.0])
        ),
        u0,
        sir_history,
        tspan,
        p;
        constant_lags = (p.τ,)
    )
    ode_f = DelayDiffEq.ODEFunctionWrapper(prob_ddae.f, prob_ddae.h)
    @test ode_f.mass_matrix == Diagonal([1.0, 1.0, 0.0])

    int = init(prob_ddae, MethodOfSteps(Rosenbrock23()))
    @test int.isdae
    @test int.f.mass_matrix == Diagonal([1.0, 1.0, 0.0])

    for (alg, reltol) in (
            (Rosenbrock23(), nothing),
            (Rodas4(), nothing),
            (QNDF(), 1.0e-6),
            (Trapezoid(), 1.0e-6),
        )
        sol_ddae = solve(prob_ddae, MethodOfSteps(alg); reltol = reltol)
        sol = appxtrue(sol_ddae, sol_dde)
        @test sol.errors[:L∞] < 5.0e-3
    end
end
