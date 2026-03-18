using StochasticDiffEq, Test
import SciMLBase
using Random

# Additive SDE test problem
f_additive(u, p, t) = @. p[2] / sqrt(1 + t) - u / (2 * (1 + t))
σ_additive(u, p, t) = @. p[1] * p[2] / sqrt(1 + t)
p = (0.1, 0.05)
additive_analytic(u0, p, t, W) = @. u0 / sqrt(1 + t) + p[2] * (t + p[1] * W) / sqrt(1 + t)
ff_additive = SDEFunction(f_additive, σ_additive, analytic = additive_analytic)
prob_sde_additive = SDEProblem(ff_additive, σ_additive, 1.0, (0.0, 1.0), p)

Random.seed!(100)

# Test default (no algorithm specified) - should use SOSRI
prob = prob_sde_additive
sol = solve(prob, dt = 1 / 2^(3))
@test sol.alg isa SOSRI

# Test with :additive hint - should use SOSRA
sol = solve(prob, dt = 1 / 2^(3), alg_hints = [:additive])
@test sol.alg isa SOSRA

# Test with :stratonovich hint - should use RKMil with Stratonovich interpretation
sol = solve(prob, dt = 1 / 2^(3), alg_hints = [:stratonovich])
@test SciMLBase.alg_interpretation(sol.alg) ==
    SciMLBase.AlgorithmInterpretation.Stratonovich
@test sol.alg isa RKMil

# Non-diagonal noise test problem
f = (du, u, p, t) -> du .= 1.01u
g = function (du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[2]
    du[2, 1] = 1.2u[1]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    return du[2, 4] = 1.8u[2]
end
prob = SDEProblem(f, g, ones(2), (0.0, 1.0), noise_rate_prototype = zeros(2, 4))

# Test default with non-diagonal noise - should use LambaEM
sol = solve(prob, dt = 1 / 2^(3))
@test sol.alg isa LambaEM

# Test with :stiff hint - should use ISSEM
sol = solve(prob, dt = 1 / 2^(3), alg_hints = [:stiff])
@test sol.alg isa ISSEM

# Test with :additive hint - should still use SOSRA (overrides non-diagonal)
sol = solve(prob, dt = 1 / 2^(3), alg_hints = [:additive])
@test sol.alg isa SOSRA

# Test with :stratonovich hint - should use LambaEulerHeun
sol = solve(prob, dt = 1 / 2^(3), alg_hints = [:stratonovich])
@test sol.alg isa LambaEulerHeun
