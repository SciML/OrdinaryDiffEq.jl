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
);

# Explicit low-order methods that don't require opts.delta
sol = solve(prob, MethodOfSteps(EM()), dt = 0.01)
@test sol.u[end] != zeros(1)
sol = solve(prob, MethodOfSteps(EulerHeun()), dt = 0.01)
@test sol.u[end] != zeros(1)

# TODO: RKMil, LambaEM, LambaEulerHeun, and other algorithms that access
# integrator.opts.delta need the full DEOptions constructor (not the legacy one
# without the delta field). This requires updating the DEOptions construction
# in the merged __init to use the non-legacy constructor for SDDE problems.
