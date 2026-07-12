using StochasticDiffEqROCK, Test, Random
import SciMLBase

# Non-diagonal (general) noise for KomBurSROCK2, see SciML/OrdinaryDiffEq.jl#3862.
# Commutative geometric SDE du_i = a_i u_i dt + Σ_k b[i,k] u_i ∘dW_k, whose exact
# Stratonovich mean is E[u_i(T)] = u0_i exp(a_i T + T/2 Σ_k b[i,k]^2).
@testset "KomBurSROCK2 non-diagonal noise" begin
    a_kb = [1.01, 1.01]
    b_kb = [0.3 0.1 0.2; 0.1 0.3 0.2]
    kb_f(u, p, t) = a_kb .* u
    kb_g(u, p, t) = b_kb .* u
    kb_f!(du, u, p, t) = (@. du = a_kb * u)
    kb_g!(du, u, p, t) = (@. du = b_kb * u)
    u0_kb = [0.5, 0.25]
    tspan_kb = (0.0, 1.0)
    prob_kb_oop = SDEProblem(kb_f, kb_g, u0_kb, tspan_kb; noise_rate_prototype = zeros(2, 3))
    prob_kb_iip = SDEProblem(kb_f!, kb_g!, u0_kb, tspan_kb; noise_rate_prototype = zeros(2, 3))

    # runs to completion with finite output (guards the DimensionMismatch / UndefVarError crashes)
    for prob_kb in (prob_kb_oop, prob_kb_iip)
        Random.seed!(100)
        sol_kb = solve(prob_kb, KomBurSROCK2(), dt = 1 / 2^6, adaptive = false)
        @test SciMLBase.successful_retcode(sol_kb)
        @test all(isfinite, sol_kb.u[end])
    end

    # weak accuracy against the exact Stratonovich mean (guards the stage sign/term math)
    T = tspan_kb[2]
    exact_mean = u0_kb .* exp.(a_kb .* T .+ (T / 2) .* vec(sum(b_kb .^ 2, dims = 2)))
    for prob_kb in (prob_kb_oop, prob_kb_iip)
        Random.seed!(200)
        N = 20000
        m = zeros(2)
        for _ in 1:N
            m .+= solve(
                prob_kb, KomBurSROCK2(); dt = 1 / 8, adaptive = false,
                save_everystep = false, save_start = false
            ).u[end]
        end
        m ./= N
        @test maximum(abs.(m .- exact_mean) ./ exact_mean) < 0.03
    end
end
