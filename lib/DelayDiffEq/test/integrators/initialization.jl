using DelayDiffEq
using OrdinaryDiffEqRosenbrock
using LinearAlgebra
using Test

@testset "CheckInit" begin
    u0_good = [0.99, 0.01, 0.0]
    sir_history(p, t) = [1.0, 0.0, 0.0]
    tspan = (0.0, 40.0)
    p = (γ = 0.5, τ = 4.0)

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
        u0_good,
        sir_history,
        tspan,
        p;
        constant_lags = (p.τ,)
    )
    alg = MethodOfSteps(Rosenbrock23())
    @test_nowarn init(prob_ddae, alg)
    prob_ddae.u0[1] = 2.0
    @test_throws SciMLBase.CheckInitFailureError init(prob_ddae, alg)
end
