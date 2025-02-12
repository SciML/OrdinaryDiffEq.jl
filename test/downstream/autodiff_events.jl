using SciMLSensitivity
using OrdinaryDiffEq, OrdinaryDiffEqCore, Calculus, Test
using Zygote

function f(du, u, p, t)
    du[1] = u[2]
    du[2] = -p[1]
end

function condition(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end

function affect!(integrator)
    @show integrator.t
    println("bounced.")
    integrator.u[2] = -integrator.p[2] * integrator.u[2]
end

cb = ContinuousCallback(condition, affect!)
p = [9.8, 0.8]
prob = ODEProblem(f, eltype(p).([1.0, 0.0]), eltype(p).((0.0, 1.0)), copy(p))

function test_f(p)
    _prob = remake(prob, p = p)
    solve(_prob, Tsit5(), abstol = 1e-14, reltol = 1e-14, callback = cb,
    save_everystep = false).u[end]
end
findiff = Calculus.finite_difference_jacobian(test_f, p)
findiff

using ForwardDiff
ad = ForwardDiff.jacobian(test_f, p)
ad

@test ad ≈ findiff

function test_f2(p, sensealg = ForwardDiffSensitivity(), controller = nothing,
        alg = Tsit5())
    _prob = remake(prob, p = p)
    u = solve(_prob, alg, sensealg = sensealg, controller = controller,
        abstol = 1e-14, reltol = 1e-14, callback = cb, save_everystep = false)
    u[end][end]
end

@test test_f2(p) == test_f(p)[end]

g1 = Zygote.gradient(θ -> test_f2(θ, ForwardDiffSensitivity()), p)
g2 = Zygote.gradient(θ -> test_f2(θ, ReverseDiffAdjoint()), p)
g3 = Zygote.gradient(θ -> test_f2(θ, ReverseDiffAdjoint(), IController()), p)
g4 = Zygote.gradient(θ -> test_f2(θ, ReverseDiffAdjoint(), PIController(7 // 50, 2 // 25)),
    p)
@test_broken g5 = Zygote.gradient(
    θ -> test_f2(θ, ReverseDiffAdjoint(),
        PIDController(1 / 18.0, 1 / 9.0, 1 / 18.0)),
    p)
g6 = Zygote.gradient(
    θ -> test_f2(θ, ForwardDiffSensitivity(),
        OrdinaryDiffEqCore.PredictiveController(), TRBDF2()),
    p)
@test_broken g7 = Zygote.gradient(
    θ -> test_f2(θ, ReverseDiffAdjoint(),
        OrdinaryDiffEqCore.PredictiveController(),
        TRBDF2()),
    p)

@test g1[1] ≈ findiff[2, 1:2]
@test g2[1] ≈ findiff[2, 1:2]
@test g3[1] ≈ findiff[2, 1:2]
@test g4[1] ≈ findiff[2, 1:2]
@test_broken g5[1] ≈ findiff[2, 1:2]
@test g6[1] ≈ findiff[2, 1:2]
@test_broken g7[1] ≈ findiff[2, 1:2]
