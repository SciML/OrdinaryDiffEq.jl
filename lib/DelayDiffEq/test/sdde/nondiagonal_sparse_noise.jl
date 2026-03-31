using DelayDiffEq
using StochasticDiffEqLowOrder
using Random
using SparseArrays
using Test

@testset "SDDE non-diagonal sparse noise" begin
    function sir_dde!(du, u, h, p, t)
        (S, I, R) = u
        (β, c, τ) = p
        N = S + I + R
        infection = β * c * I / N * S
        (Sd, Id, Rd) = h(p, t - τ)
        Nd = Sd + Id + Rd
        recovery = β * c * Id / Nd * Sd
        @inbounds begin
            du[1] = -infection
            du[2] = infection - recovery
            du[3] = recovery
        end
        return nothing
    end

    A = zeros(3, 2)
    A[1, 1] = 1
    A[2, 1] = 1
    A[2, 2] = 1
    A[3, 2] = 1
    A = SparseArrays.sparse(A)

    function sir_delayed_noise!(du, u, h, p, t)
        (S, I, R) = u
        (β, c, τ) = p
        N = S + I + R
        infection = β * c * I / N * S
        (Sd, Id, Rd) = h(p, t - τ)
        Nd = Sd + Id + Rd
        recovery = β * c * Id / Nd * Sd
        du[1, 1] = -sqrt(max(0.0, infection))
        du[2, 1] = sqrt(max(0.0, infection))
        du[2, 2] = -sqrt(max(0.0, recovery))
        return du[3, 2] = sqrt(max(0.0, recovery))
    end

    function condition(u, t, integrator)
        return u[2]
    end
    function affect!(integrator)
        return integrator.u[2] = 0.0
    end
    cb = ContinuousCallback(condition, affect!)

    δt = 0.1
    tmax = 40.0
    tspan = (0.0, tmax)
    u0 = [990.0, 10.0, 0.0]

    function sir_history(p, t)
        return [1000.0, 0.0, 0.0]
    end

    p = [0.05, 10.0, 4.0]
    Random.seed!(1234)

    prob_sdde = SDDEProblem(
        sir_dde!, sir_delayed_noise!, u0, sir_history, tspan, p;
        noise_rate_prototype = A
    )
    # Use EM instead of LambaEM (LambaEM needs opts.delta support)
    sol_sdde = solve(prob_sdde, MethodOfSteps(EM()), dt = 0.1, callback = cb)
    @test sol_sdde.retcode == ReturnCode.Success
    @test length(sol_sdde.t) > 10
end
