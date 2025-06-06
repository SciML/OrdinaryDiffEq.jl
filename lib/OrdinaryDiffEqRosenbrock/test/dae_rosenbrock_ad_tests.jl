using OrdinaryDiffEqRosenbrock, LinearAlgebra, ForwardDiff, Test
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit, ShampineCollocationInit
using ADTypes: AutoForwardDiff, AutoFiniteDiff

afd_cs3 = AutoForwardDiff(chunksize=3)
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
function rober(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    [-k₁ * y₁ + k₃ * y₂ * y₃,
        k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2,
        y₁ + y₂ + y₃ - 1]
end
M = Diagonal([1.0, 1.0, 0.0])
roberf = ODEFunction{true, SciMLBase.AutoSpecialize}(rober, mass_matrix = M)
roberf_oop = ODEFunction{false, SciMLBase.AutoSpecialize}(rober, mass_matrix = M)
prob_mm = ODEProblem(roberf, [1.0, 0.0, 0.2], (0.0, 1e5), (0.04, 3e7, 1e4))
prob_mm_oop = ODEProblem(roberf_oop, [1.0, 0.0, 0.2], (0.0, 1e5), (0.04, 3e7, 1e4))
# Both should be inferrable so long as AutoSpecialize is used...
@test_broken sol = @inferred solve(prob_mm, Rodas5P(), reltol = 1e-8, abstol = 1e-8) 
sol = @inferred solve(prob_mm_oop, Rodas5P(), reltol = 1e-8, abstol = 1e-8)

# These tests flex differentiation of the solver and through the initialization
# To only test the solver part and isolate potential issues, set the initialization to consistent
@testset "Inplace: $(isinplace(_prob)), BrownBasic: $(initalg isa BrownFullBasicInit), Autodiff: $autodiff" for _prob in [
        prob_mm, prob_mm_oop],
    initalg in [BrownFullBasicInit(), ShampineCollocationInit()], autodiff in [AutoForwardDiff(chunksize=3), AutoFiniteDiff()]

    alg = Rodas5P(; autodiff)
    function f(p)
        sol = @inferred solve(remake(_prob, p = p), alg, abstol = 1e-14,
            reltol = 1e-14, initializealg = initalg)
        sum(sol)
    end
    @test ForwardDiff.gradient(f, [0.04, 3e7, 1e4])≈[0, 0, 0] atol=1e-8
end
