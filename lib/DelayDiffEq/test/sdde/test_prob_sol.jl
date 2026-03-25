using DelayDiffEq
using StochasticDiffEqLowOrder
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
padd = [1.0, -4.0, -2.0, 10.0, -0.0, -0.0, 0.1]

prob = SDDEProblem(
    hayes_modelf, hayes_modelg, [1.0], h, tspan, pmul;
    constant_lags = (pmul[1],)
)

# Explicit low-order methods
sol = solve(prob, MethodOfSteps(EM()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(LambaEM()), dt = 0.01, adaptive = false)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(EulerHeun()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(LambaEulerHeun()), dt = 0.01, adaptive = false)
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

# TODO: Additive noise methods (SRA, SRA1, SOSRA, SOSRA2) need
# alg_needs_extra_process support in _create_sdde_noise for the dZ process
