using DelayDiffEq
using StochasticDiffEqLowOrder
using StochasticDiffEqHighOrder
using Test

# Hayes Equation
begin
    function hayes_modelf(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        return du .= a .* u .+ b .* h(p, t - τ) .+ c
    end
    function hayes_modelg(du, u, h, p, t)
        τ, a, b, c, α, β, γ = p
        return du .= α .* u .+ β .* h(p, t - τ) .+ γ
    end
    h(p, t) = (ones(1) .+ t)
    tspan = (0.0, 0.1)
end

pmul = [1.0, -4.0, -2.0, 10.0, -1.3, -1.2, 1.1]
p_additive = [1.0, -4.0, -2.0, 10.0, -0.0, -0.0, 0.1]

prob = SDDEProblem(
    hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
    constant_lags = (pmul[1],)
)

# Explicit low-order methods
sol = solve(prob, MethodOfSteps(EM()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(LambaEM()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(EulerHeun()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(LambaEulerHeun()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(RKMil()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(
    prob, MethodOfSteps(RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich)),
    dt = 0.01, adaptive = false
)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(RKMilCommute()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(
    prob, MethodOfSteps(RKMilCommute(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich)),
    dt = 0.01, adaptive = false
)
@test sol.u[end] != zeros(1)

# High-order SRI methods (need dZ extra process)
sol = solve(prob, MethodOfSteps(SRI()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(SRIW1()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(SOSRI()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(SOSRI2()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)

# Additive noise problems (need dZ extra process)
prob_add = SDDEProblem(
    hayes_modelf, hayes_modelg, [1.0], h, tspan, p_additive;
    constant_lags = (p_additive[1],)
)
sol = solve(prob_add, MethodOfSteps(SRA()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob_add, MethodOfSteps(SRA1()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob_add, MethodOfSteps(SOSRA()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob_add, MethodOfSteps(SOSRA2()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
