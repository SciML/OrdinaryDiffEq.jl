# Version-dependent AD backend selection
# Enzyme: all AD test lanes
# Zygote: Julia <= 1.11 only
# Mooncake: works on all Julia versions for the ForwardDiffSensitivity-based
#   gradient tests below via SciMLSensitivityMooncakeExt. Mooncake does NOT
#   support ReverseDiffAdjoint (TypeError on the ODESolution CoDual), so the
#   ReverseDiffAdjoint-based gradient tests are not run with Mooncake.
# ForwardDiff: all versions

const JULIA_VERSION_ALLOWS_ZYGOTE = VERSION < v"1.12" && isempty(VERSION.prerelease)

using SciMLSensitivity
using OrdinaryDiffEq, OrdinaryDiffEqCore, Test
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore: IController, PIController, PIDController
using ADTypes
import DifferentiationInterface as DI
import ReverseDiff
using Enzyme
using Mooncake  # Load Mooncake after DI to ensure extension is loaded

# Load version-dependent packages
if JULIA_VERSION_ALLOWS_ZYGOTE
    using Zygote
    get_gradient_backends() = [AutoZygote()]
else
    # On Julia 1.12+, Zygote isn't available; the gradient_backends list
    # is empty here because the existing gradient testset exercises sensealgs
    # (ReverseDiffAdjoint, PIDController, ...) that Mooncake does not support.
    # Mooncake-specific gradient tests for the sensealgs Mooncake DOES support
    # are added in a separate testset below.
    get_gradient_backends() = []
end
get_jacobian_backends() = [AutoForwardDiff()]

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
# There is exactly one bounce on this interval, so the final-state Jacobian has
# a closed form that avoids finite-difference noise from perturbing the root time.
impact_speed = sqrt(2 * p[1])
expected = [
    (1 + p[2]) / impact_speed - 1 / 2 impact_speed - 2
    (1 + p[2]) / impact_speed - 1 impact_speed
]

# Test jacobians with all available backends
@testset "Jacobian tests" begin
    for backend in get_jacobian_backends()
        @testset "$(typeof(backend))" begin
            ad = DI.jacobian(test_f, backend, p)
            @test ad ≈ expected
        end
    end
end

# Enzyme fails on ContinuousCallback with "mixed activity for jl_new_struct"
@testset "Enzyme callback limitation (jacobian)" begin
    @test_broken (
        ad = DI.jacobian(test_f, AutoEnzyme(mode = Enzyme.set_runtime_activity(Enzyme.Forward)), p);
        ad ≈ expected
    )
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
    return u[end]
end

@test test_f2(p) == test_f(p)[end]

@testset "Out-of-place Rosenbrock23 with ReverseDiffAdjoint" begin
    t_save = 0.0:0.05:0.1
    oop_prob = ODEProblem{false}(
        (u, p, t) -> [u[2], -p[1]], [1.0, 0.0], (0.0, 0.1), [9.81, 1.0]
    )
    alg = Rosenbrock23(autodiff = AutoFiniteDiff())

    function oop_loss(p)
        sol = solve(
            oop_prob; p, alg, saveat = t_save, sensealg = ReverseDiffAdjoint(),
            abstol = 1.0e-10, reltol = 1.0e-10
        )
        return sum(sol[1, :])
    end

    @test ReverseDiff.gradient(oop_loss, [9.81, 1.0]) ≈ [-0.00625, 0.0]
end

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
    g5 = DI.gradient(
        θ -> test_f2(
            θ, ReverseDiffAdjoint(),
            PIDController(1 / 18.0, 1 / 9.0, 1 / 18.0)
        ),
        backend,
        p
    )
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

    @test g1 ≈ expected[2, 1:2]
    @test g2 ≈ expected[2, 1:2]
    @test g3 ≈ expected[2, 1:2]
    @test g4 ≈ expected[2, 1:2]
    @test g5 ≈ expected[2, 1:2]
    @test g6 ≈ expected[2, 1:2]
    @test_broken g7 ≈ expected[2, 1:2]
end

# Enzyme fails on ContinuousCallback with "mixed activity for jl_new_struct"
@testset "Enzyme callback limitation (gradient)" begin
    @test (
        g = DI.gradient(θ -> test_f2(θ, ForwardDiffSensitivity()), AutoEnzyme(mode = Enzyme.set_runtime_activity(Enzyme.Reverse)), p);
        g ≈ expected[2, 1:2]
    )
end

# Mooncake gradient tests (all Julia versions). Only the sensealgs Mooncake
# actually supports are exercised here:
#   - ForwardDiffSensitivity: works for the standard, IController, and
#     PredictiveController/TRBDF2 cases.
#   - ReverseDiffAdjoint, PI/PIDController-driven solves: NOT exercised because
#     Mooncake currently throws a TypeError on the resulting ODESolution CoDual.
#
# `test_f2` returns `sol[end][end]`. The `sol[Int]` rrule in
# `SciMLBaseMooncakeExt._scatter_pullback` is currently broken: it treats the
# integer as a variable index, but for `ODESolution` `sol[Int]` returns the
# state at that time index, leading to a BoundsError. As a workaround we
# define a Mooncake-only variant that reaches into `sol.u` directly to avoid
# the SciMLBase getindex rrule entirely.
function test_f2_mc(p, sensealg, controller = nothing, alg = Tsit5())
    _prob = remake(prob, p = p)
    sol = solve(
        _prob, alg, sensealg = sensealg, controller = controller,
        abstol = 1.0e-14, reltol = 1.0e-14, callback = cb, save_everystep = false
    )
    return sol.u[end][end]
end

@testset "Mooncake gradient tests" begin
    backend = AutoMooncake(; config = nothing)
    g1 = DI.gradient(θ -> test_f2_mc(θ, ForwardDiffSensitivity()), backend, p)
    g6 = DI.gradient(
        θ -> test_f2_mc(
            θ, ForwardDiffSensitivity(),
            OrdinaryDiffEqCore.PredictiveController(), TRBDF2()
        ),
        backend,
        p
    )
    @test g1 ≈ expected[2, 1:2] rtol = 1.0e-5
    @test g6 ≈ expected[2, 1:2] rtol = 1.0e-5
end
