# Version-dependent AD backend selection
# Enzyme/Zygote: Julia <= 1.11 only (see https://github.com/EnzymeAD/Enzyme.jl/issues/2699)
# Mooncake: gradient support for ODE solves is currently broken (see discrete_adjoints.jl)
# ForwardDiff: all versions

const JULIA_VERSION_ALLOWS_ENZYME_ZYGOTE = VERSION < v"1.12" && isempty(VERSION.prerelease)

using SciMLSensitivity
using OrdinaryDiffEq, OrdinaryDiffEqCore, FiniteDiff, Test
using ADTypes
import DifferentiationInterface as DI

# Load version-dependent packages
if JULIA_VERSION_ALLOWS_ENZYME_ZYGOTE
    using Enzyme
    using Zygote
    get_gradient_backends() = [AutoZygote()]
    get_jacobian_backends() = [AutoForwardDiff()]
else
    # On Julia 1.12+, skip gradient tests since Zygote/Enzyme aren't available
    # and Mooncake gradient support for ODE solves is broken
    get_gradient_backends() = []
    get_jacobian_backends() = [AutoForwardDiff()]
end

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

# Test jacobians with all available backends
@testset "Jacobian tests" begin
    for backend in get_jacobian_backends()
        @testset "$(typeof(backend))" begin
            ad = DI.jacobian(test_f, backend, p)
            @test ad ≈ findiff
        end
    end
end

# Enzyme fails on ContinuousCallback with "mixed activity for jl_new_struct"
if JULIA_VERSION_ALLOWS_ENZYME_ZYGOTE
    @testset "Enzyme callback limitation (jacobian)" begin
        @test_broken (
            ad = DI.jacobian(test_f, AutoEnzyme(mode = Enzyme.set_runtime_activity(Enzyme.Forward)), p);
            ad ≈ findiff
        )
    end
end

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

# Test gradients with all available reverse-mode backends
@testset "Gradient tests with $backend" for backend in get_gradient_backends()
    g1 = DI.gradient(θ -> test_f2(θ, ForwardDiffSensitivity()), backend, p)
    g2 = DI.gradient(θ -> test_f2(θ, ReverseDiffAdjoint()), backend, p)
    g3 = DI.gradient(θ -> test_f2(θ, ReverseDiffAdjoint(), IController()), backend, p)
    g4 = DI.gradient(
        θ -> test_f2(θ, ReverseDiffAdjoint(), PIController(7 // 50, 2 // 25)),
        backend,
        p
    )
    if VERSION >= v"1.11"
        g5 = DI.gradient(
            θ -> test_f2(
                θ, ReverseDiffAdjoint(),
                PIDController(1 / 18.0, 1 / 9.0, 1 / 18.0)
            ),
            backend,
            p
        )
    else
        @test_broken g5 = DI.gradient(
            θ -> test_f2(
                θ, ReverseDiffAdjoint(),
                PIDController(1 / 18.0, 1 / 9.0, 1 / 18.0)
            ),
            backend,
            p
        )
    end
    g6 = DI.gradient(
        θ -> test_f2(
            θ, ForwardDiffSensitivity(),
            OrdinaryDiffEqCore.PredictiveController(), TRBDF2()
        ),
        backend,
        p
    )
    @test_broken g7 = DI.gradient(
        θ -> test_f2(
            θ, ReverseDiffAdjoint(),
            OrdinaryDiffEqCore.PredictiveController(),
            TRBDF2()
        ),
        backend,
        p
    )

    @test g1 ≈ findiff[2, 1:2]
    @test g2 ≈ findiff[2, 1:2]
    @test g3 ≈ findiff[2, 1:2]
    @test g4 ≈ findiff[2, 1:2]
    if VERSION >= v"1.11"
        @test g5 ≈ findiff[2, 1:2]
    end
    @test g6 ≈ findiff[2, 1:2]
    @test_broken g7 ≈ findiff[2, 1:2]
end

# Enzyme fails on ContinuousCallback with "mixed activity for jl_new_struct"
if JULIA_VERSION_ALLOWS_ENZYME_ZYGOTE
    @testset "Enzyme callback limitation (gradient)" begin
        @test_broken (
            g = DI.gradient(θ -> test_f2(θ, ForwardDiffSensitivity()), AutoEnzyme(mode = Enzyme.set_runtime_activity(Enzyme.Reverse)), p);
            g ≈ findiff[2, 1:2]
        )
    end
end
