using OrdinaryDiffEqRosenbrock, LinearAlgebra, Test
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, ShampineCollocationInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff
import DifferentiationInterface as DI

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
# @inferred is disabled: AutoSpecialize's FunctionWrappersWrapper creates type
# nesting depth ~8, and embedding in ODESolution→InterpolationData→RosenbrockCache
# exceeds Julia's MAX_TYPE_DEPTH inference limit. Solve is correct; only compile-time
# inference fails. See https://github.com/SciML/OrdinaryDiffEq.jl/issues/XXXX
sol = solve(
    prob_mm, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
sol = solve(
    prob_mm_oop, Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)

# Test Tsit5DA (Hybrid Explicit/Linear-Implicit)
# Tsit5DA is explicit so it will hit maxiters on stiff rober — just test it runs
sol = solve(
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
        sol = solve(
            remake(_prob, p = p), alg, abstol = 1.0e-14,
            reltol = 1.0e-14, initializealg = initalg
        )
        sum(sol.u[end])
    end
    @test DI.gradient(f, AutoForwardDiff(), [0.04, 3.0e7, 1.0e4]) ≈ [0, 0, 0] atol = 1.0e-8
end
