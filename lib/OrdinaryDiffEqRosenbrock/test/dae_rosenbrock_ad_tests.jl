using OrdinaryDiffEqRosenbrock, LinearAlgebra, Test
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, ShampineCollocationInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
import DifferentiationInterface as DI
using SciMLBase: ODEProblem
using ForwardDiff

afd_cs3 = AutoForwardDiff(chunksize = 3)
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    return [
        -k₁ * y₁ + k₃ * y₂ * y₃,
        k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
        y₁ + y₂ + y₃ - 1,
    ]
end
M = Diagonal([1.0, 1.0, 0.0])
roberf = ODEFunction{true, SciMLBase.AutoSpecialize}(rober, mass_matrix = M)
roberf_oop = ODEFunction{false, SciMLBase.AutoSpecialize}(rober, mass_matrix = M)
prob_mm = ODEProblem(roberf, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
prob_mm_oop = ODEProblem(roberf_oop, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
# Both should be inferable so long as AutoSpecialize is used...
sol = @inferred solve(
    prob_mm, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
sol = @inferred solve(
    prob_mm_oop, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)

# Test Tsit5DA (Hybrid Explicit/Linear-Implicit)
# Tsit5DA is explicit so it will hit maxiters on stiff rober — just test inference
sol = @inferred solve(
    prob_mm, Tsit5DA(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)

# These tests flex differentiation of the solver and through the initialization
# To only test the solver part and isolate potential issues, set the initialization to consistent
@testset "IIP=$(isinplace(_prob)) Brown=$(initalg isa BrownFullBasicInit) AD=$autodiff Alg=$AlgName" for _prob in [
            prob_mm, prob_mm_oop,
        ],
        initalg in [BrownFullBasicInit(), ShampineCollocationInit()],
        autodiff in [AutoForwardDiff(chunksize = 3), AutoFiniteDiff()],
        (Alg, AlgName) in [(Rodas5P, "Rodas5P"), (Tsit5DA, "Tsit5DA")]

    alg = Alg(; autodiff)
    function f(p)
        sol = @inferred solve(
            remake(_prob, p = p), alg, abstol = 1.0e-14,
            reltol = 1.0e-14, initializealg = initalg
        )
        sum(sol.u[end])
    end
    @test DI.gradient(f, AutoForwardDiff(), [0.04, 3.0e7, 1.0e4]) ≈ [0, 0, 0] atol = 1.0e-8
end

# Regression test for issue #3486 — nested ForwardDiff (hessian / gradient-of-gradient)
# through a Rosenbrock solve used to error with:
#   MethodError: no method matching Float64(::ForwardDiff.Dual{...})
# ...when the Jacobian-reuse state tried to unwrap one Dual level via
# `ForwardDiff.value(dtgamma)` and assign it to a strictly-typed `::Float64`
# field. Unconditional unwrapping collapses the still-active inner derivative,
# so the fix stores `dtgamma` verbatim (Dual and all) in `JacReuseState`.
@testset "Second-order ForwardDiff through solve (issue #3486)" begin
    rhs_oop_simple(u, p, t) = [-p[1] * u[1] * u[1]]
    function loss_oop(p)
        prob = ODEProblem{false}(rhs_oop_simple, [1.0], (0.0, p[2]), p)
        sol = solve(prob, Rodas5P())
        return only(last(sol.u))
    end
    g = ForwardDiff.gradient(loss_oop, [1.0, 1.0])
    @test all(isfinite, g)
    H = ForwardDiff.hessian(loss_oop, [1.0, 1.0])
    @test all(isfinite, H)
end
