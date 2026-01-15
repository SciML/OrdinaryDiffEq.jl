# Skip Enzyme tests on Julia 1.12+ prerelease versions
if !isempty(VERSION.prerelease)
    @warn "Skipping Enzyme tests on Julia prerelease version $(VERSION)"
    exit(0)
end

# Enzyme and Zygote are only supported on Julia <= 1.11
if VERSION >= v"1.12"
    @warn "Skipping Enzyme/Zygote tests on Julia $(VERSION) - only supported on Julia <= 1.11"
    exit(0)
end

using SciMLSensitivity
using OrdinaryDiffEq, OrdinaryDiffEqCore, FiniteDiff, Test
using ADTypes
import DifferentiationInterface as DI
using Enzyme
using Zygote

function f(du, u, p, t)
    du[1] = u[2]
    return du[2] = -p[1]
end

function condition(u, t, integrator) # Event when event_f(u,t) == 0
    return u[1]
end

function affect!(integrator)
    @show integrator.t
    println("bounced.")
    return integrator.u[2] = -integrator.p[2] * integrator.u[2]
end

cb = ContinuousCallback(condition, affect!)
p = [9.8, 0.8]
prob = ODEProblem(f, eltype(p).([1.0, 0.0]), eltype(p).((0.0, 1.0)), copy(p))

function test_f(p)
    _prob = remake(prob, p = p)
    return solve(
        _prob, Tsit5(), abstol = 1.0e-14, reltol = 1.0e-14, callback = cb,
        save_everystep = false
    ).u[end]
end
findiff = FiniteDiff.finite_difference_jacobian(test_f, p)
findiff

# Test with DifferentiationInterface using AutoForwardDiff
ad = DI.jacobian(test_f, AutoForwardDiff(), p)
ad

@test ad ≈ findiff

# Test with Enzyme forward mode
ad_enzyme = DI.jacobian(test_f, AutoEnzyme(mode = Enzyme.Forward), p)
@test ad_enzyme ≈ findiff

function test_f2(
        p, sensealg = ForwardDiffSensitivity(), controller = nothing,
        alg = Tsit5()
    )
    _prob = remake(prob, p = p)
    u = solve(
        _prob, alg, sensealg = sensealg, controller = controller,
        abstol = 1.0e-14, reltol = 1.0e-14, callback = cb, save_everystep = false
    )
    return u[end][end]
end

@test test_f2(p) == test_f(p)[end]

# Use DifferentiationInterface with AutoZygote backend for gradient computations
g1 = DI.gradient(θ -> test_f2(θ, ForwardDiffSensitivity()), AutoZygote(), p)
g2 = DI.gradient(θ -> test_f2(θ, ReverseDiffAdjoint()), AutoZygote(), p)
g3 = DI.gradient(θ -> test_f2(θ, ReverseDiffAdjoint(), IController()), AutoZygote(), p)
g4 = DI.gradient(
    θ -> test_f2(θ, ReverseDiffAdjoint(), PIController(7 // 50, 2 // 25)),
    AutoZygote(),
    p
)
@test_broken g5 = DI.gradient(
    θ -> test_f2(
        θ, ReverseDiffAdjoint(),
        PIDController(1 / 18.0, 1 / 9.0, 1 / 18.0)
    ),
    AutoZygote(),
    p
)
g6 = DI.gradient(
    θ -> test_f2(
        θ, ForwardDiffSensitivity(),
        OrdinaryDiffEqCore.PredictiveController(), TRBDF2()
    ),
    AutoZygote(),
    p
)
@test_broken g7 = DI.gradient(
    θ -> test_f2(
        θ, ReverseDiffAdjoint(),
        OrdinaryDiffEqCore.PredictiveController(),
        TRBDF2()
    ),
    AutoZygote(),
    p
)

@test g1 ≈ findiff[2, 1:2]
@test g2 ≈ findiff[2, 1:2]
@test g3 ≈ findiff[2, 1:2]
@test g4 ≈ findiff[2, 1:2]
@test_broken g5 ≈ findiff[2, 1:2]
@test g6 ≈ findiff[2, 1:2]
@test_broken g7 ≈ findiff[2, 1:2]
