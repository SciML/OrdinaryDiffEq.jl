using StochasticDiffEq, OrdinaryDiffEq, Test
import SciMLBase

function f(du, u, p, t)
    return du[1] = u[1]
end
function g(du, u, p, t)
    return 0.0
end

u0 = [1.0]
prob = SDEProblem{true}(f, g, u0, (0.0, 0.1))

sol_ito = solve(prob, RKMil{SciMLBase.AlgorithmInterpretation.Ito}())
@test length(sol_ito) < 100

sol_strato = solve(prob, RKMil{SciMLBase.AlgorithmInterpretation.Stratonovich}(); dt = 1.0e-2)
@test length(sol_strato) < 100

sol_ito = solve(prob, RKMil())
@test length(sol_ito) < 100

sol_strato = solve(prob, RKMil(interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich); dt = 1.0e-2)
@test length(sol_strato) < 100

sol_leh = solve(prob, LambaEulerHeun())
@test length(sol_leh) < 100

# Test that SDE solvers with zero noise converge to the same solution as ODE solvers
# Regression test for issue #636: ImplicitRKMil error estimation bug
# https://github.com/SciML/StochasticDiffEq.jl/issues/636

# Linear decay system: du/dt = -λu with λ = 1 (mildly stiff)
# Analytic solution: u(t) = u0 * exp(-λt)
function linear_decay!(du, u, p, t)
    λ = p
    du[1] = -λ * u[1]
    return du[2] = -λ * u[2]
end

function zero_noise_linear!(du, u, p, t)
    return du .= 0
end

# Parameters - use λ=1 for mild stiffness (faster tests)
λ = 1.0
p_decay = λ
u0_decay = [1.0, 2.0]
tspan_decay = (0.0, 1.0)

# Create SDE and ODE problems
sde_prob_decay = SDEProblem(
    linear_decay!, zero_noise_linear!, u0_decay, tspan_decay, p_decay
)
ode_prob_decay = ODEProblem(linear_decay!, u0_decay, tspan_decay, p_decay)

# Analytic solution for reference
analytic_endpoint = u0_decay .* exp(-λ * tspan_decay[2])

# Solve with ODE solver to verify setup
abstol_decay = 1.0e-4
reltol_decay = 1.0e-4
sol_ode_decay = solve(ode_prob_decay, Rodas5(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_ode_decay.u[end], analytic_endpoint, rtol = 1.0e-3)

# Test implicit SDE solvers that were fixed in this PR
# These should now converge to the same solution as the ODE solver

# ImplicitRKMil (the main solver from issue #636)
sol_implicitRKMil = solve(sde_prob_decay, ImplicitRKMil(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_implicitRKMil.u[end], analytic_endpoint, rtol = 0.05)

# ImplicitEM
sol_implicitEM = solve(sde_prob_decay, ImplicitEM(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_implicitEM.u[end], analytic_endpoint, rtol = 0.05)

# ImplicitEulerHeun
sol_implicitEulerHeun = solve(sde_prob_decay, ImplicitEulerHeun(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_implicitEulerHeun.u[end], analytic_endpoint, rtol = 0.05)

# Test other SDE solvers for consistency

# SKenCarp (stiff SDE solver)
sol_SKenCarp = solve(sde_prob_decay, SKenCarp(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_SKenCarp.u[end], analytic_endpoint, rtol = 0.05)

# ISSEM (semi-implicit)
sol_ISSEM = solve(sde_prob_decay, ISSEM(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_ISSEM.u[end], analytic_endpoint, rtol = 0.05)

# SOSRI (default adaptive SDE solver, known to work correctly)
sol_SOSRI = solve(sde_prob_decay, SOSRI(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_SOSRI.u[end], analytic_endpoint, rtol = 0.05)

# SOSRA
sol_SOSRA = solve(sde_prob_decay, SOSRA(), abstol = abstol_decay, reltol = reltol_decay)
@test isapprox(sol_SOSRA.u[end], analytic_endpoint, rtol = 0.05)
