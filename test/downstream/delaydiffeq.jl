using DelayDiffEq, DDEProblemLibrary, ADTypes
using Test

@testset "Constant delays" begin
    prob = DDEProblemLibrary.prob_dde_constant_2delays_ip
    prob_scalar = DDEProblemLibrary.prob_dde_constant_2delays_scalar
    # disable reuse
    nlsolve = nlsolve = NLNewton(fast_convergence_cutoff = 0)

    algdict = Dict(
        BS3() => 2.4e-6,
        Tsit5() => 1.0e-2,
        RK4() => 2.0e-4,
        Vern6() => 3.0e-3,
        SDIRK2(nlsolve = nlsolve) => 3.8e-1,
        TRBDF2(nlsolve = nlsolve) => 6.2e-2,
        KenCarp4(nlsolve = nlsolve) => 7.3e-2,
        Rosenbrock23() => 6.5e-4,
        Rodas4() => 2.0e-3
    )

    for (alg, error) in algdict
        ddealg = MethodOfSteps(alg)

        sol = solve(prob, ddealg)
        @test sol.errors[:l∞] < error

        sol_scalar = solve(prob_scalar, ddealg)
        @test sol.t ≈ sol_scalar.t atol = 1.0e-3
        @test sol[1, :] ≈ sol_scalar.u atol = 1.0e-3
    end
end

function lotka_volterra!(du, u, h, p, t)
    🐰, 🐺 = u
    α, β, γ, δ, τ = p
    🕥🐰 = h(p, t - τ; idxs = 1)
    du[1] = d🐰 = α * 🕥🐰 - β * 🐺 * 🐰
    du[2] = d🐺 = γ * 🐺 * 🐰 - δ * 🐺
    return nothing
end
uₒ = [1.0, 1.0]
tspan = (0.0, 10.0)
h(p, t) = [1.0, 1.0]
h(p, t; idxs = 1) = 1.0
p = [1.5, 1.0, 3.0, 1.0, 1.0]
prob = DDEProblem(lotka_volterra!, uₒ, h, tspan, p, constant_lags = (p[end],))
sol = solve(prob, MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff = AutoFiniteDiff()))))
