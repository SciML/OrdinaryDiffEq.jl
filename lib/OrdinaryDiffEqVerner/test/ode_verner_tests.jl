using OrdinaryDiffEqVerner, OrdinaryDiffEqCore, OrdinaryDiffEqExplicitRK
using DiffEqDevTools, Test
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear

# Problem mappings
const probbig = prob_ode_bigfloat2Dlinear
const probnumbig = prob_ode_bigfloatlinear
const probnum = prob_ode_linear
const prob = prob_ode_2Dlinear

testTol = 2.0

# Custom ODE functions for testing
function f!(du, u, p, t)
    return du[1] = -u[1]
end

function f(u, p, t)
    return -u
end

t_end = 64.0
setprecision(256)
prob_oop = ODEProblem(ODEFunction(f; analytic = (u0, p, t) -> exp(-t)), 1.0, (0.0, t_end))
prob_iip = ODEProblem(ODEFunction(f!; analytic = (u0, p, t) -> [exp(-t)]), [1.0], (0.0, t_end))

# -------------------------------------------------------------
# Helper for convergence test consistency
function check_convergence(dts, prob, alg, order_expected)
    # Use DiffEqDevTools for overall convergence estimate
    sim = test_convergence(dts, prob, alg)
    return @test (sim.𝒪est[:final] > order_expected) || (abs(sim.𝒪est[:final] - order_expected) < testTol)

end
# -------------------------------------------------------------

### Vern6
println("Vern6")
dts = (1 / 2) .^ (8:-1:5)
check_convergence(dts, probnumbig, Vern6(), 6)
check_convergence(dts, probbig, Vern6(), 6)

tabalg = ExplicitRK(tableau = constructVernerEfficient6(BigFloat))
sol1 = solve(probnumbig, Vern6(); dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg; dt = 1 / 2^6, adaptive = false, save_everystep = false)
@test sol1.u[end] - sol2.u[end] < 1.0e-10

sol1 = solve(probbig, Vern6(); dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg; dt = 1 / 2^3, adaptive = false, save_everystep = false)
@test minimum(sol1.u[end] - sol2.u[end] .< 1.0e-10)

sol1 = solve(probbig, tabalg; dt = 1 / 2^6)
sol2 = solve(probbig, Vern6(); dt = 1 / 2^6)
@test length(sol1) == length(sol2)
@test SciMLBase.successful_retcode(sol1)
@test SciMLBase.successful_retcode(sol2)

# -------------------------------------------------------------
### Vern7
println("Vern7")
dts = (1 / 2) .^ (6:-1:3)
check_convergence(dts, probnumbig, Vern7(), 7)
check_convergence(dts, probbig, Vern7(), 7)

tabalg = ExplicitRK(tableau = constructVerner7(BigFloat))
sol1 = solve(probnumbig, Vern7(); dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg; dt = 1 / 2^6, adaptive = false, save_everystep = false)
@test sol1.u[end] - sol2.u[end] < 1.0e-10

sol1 = solve(probbig, Vern7(); dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg; dt = 1 / 2^3, adaptive = false, save_everystep = false)
@test minimum(sol1.u[end] - sol2.u[end] .< 1.0e-10)

sol1 = solve(probbig, tabalg; dt = 1 / 2^6)
sol2 = solve(probbig, Vern7(); dt = 1 / 2^6)
@test length(sol1) == length(sol2)
@test SciMLBase.successful_retcode(sol1)
@test SciMLBase.successful_retcode(sol2)

# -------------------------------------------------------------
### Vern8
println("Vern8")
dts = (1 / 2) .^ (6:-1:3)
check_convergence(dts, probnumbig, Vern8(), 8)
check_convergence(dts, probbig, Vern8(), 8)

tabalg = ExplicitRK(tableau = constructVerner8(BigFloat))
sol1 = solve(probnumbig, Vern8(); dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg; dt = 1 / 2^6, adaptive = false, save_everystep = false)
@test sol1.u[end] - sol2.u[end] < 1.0e-10

sol1 = solve(probbig, Vern8(); dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg; dt = 1 / 2^3, adaptive = false, save_everystep = false)
@test minimum(sol1.u[end] - sol2.u[end] .< 1.0e-10)

sol1 = solve(prob, tabalg; dt = 1 / 2^6)
sol2 = solve(prob, Vern8(); dt = 1 / 2^6)
@test length(sol1) == length(sol2)
@test SciMLBase.successful_retcode(sol1)
@test SciMLBase.successful_retcode(sol2)

# -------------------------------------------------------------
### Vern9
println("Vern9")
dts = (1 / 2) .^ (6:-1:3)
check_convergence(dts, probnumbig, Vern9(), 9)
check_convergence(dts, probbig, Vern9(), 9)

tabalg = ExplicitRK(tableau = constructVernerEfficient9(BigFloat))
sol1 = solve(probnumbig, Vern9(); dt = 1 / 2^6, adaptive = false, save_everystep = false)
sol2 = solve(probnumbig, tabalg; dt = 1 / 2^6, adaptive = false, save_everystep = false)
@test abs(sol1.u[end] - sol2.u[end]) < 1.0e-15

sol1 = solve(probbig, Vern9(); dt = 1 / 2^3, adaptive = false, save_everystep = false)
sol2 = solve(probbig, tabalg; dt = 1 / 2^3, adaptive = false, save_everystep = false)
@test minimum(abs.(sol1.u[end] - sol2.u[end]) .< 1.0e-15)

sol1 = solve(probbig, tabalg; dt = 1 / 2^6)
sol2 = solve(probbig, Vern9(); dt = 1 / 2^6)
@test length(sol1) == length(sol2)
@test SciMLBase.successful_retcode(sol1)
@test SciMLBase.successful_retcode(sol2)

# -------------------------------------------------------------
### RKV76IIa
println("RKV76IIa")
# dts = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 0.0009765625, 0.00048828125, 0.000244140625, 0.0001220703125, 0.00006103515625]
dts = [2, 1, 0.5, 0.25]

check_convergence(dts, prob_oop, RKV76IIa(), 7)
