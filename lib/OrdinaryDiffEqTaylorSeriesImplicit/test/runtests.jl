using OrdinaryDiffEqTaylorSeriesImplicit, ODEProblemLibrary, DiffEqDevTools
using OrdinaryDiffEqRosenbrock
using Test, Random

testTol = 0.2

@testset "Implicit Solver Convergence Tests" begin
    prob_stiff = ODEProblemLibrary.prob_ode_hires
    dts = 2.0 .^ (-8:-5)
    ref_setup = Dict(:alg => Rodas5P(), :reltol => 1.0e-16, :abstol => 1.0e-16)

    # Taylor-Gauss
    sim1v1 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(1), order_q = Val(1)), ref_setup)
    @test sim1v1.𝒪est[:final] ≈ 2 atol = testTol

    sim2v2 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(2), order_q = Val(2)), ref_setup)
    @test sim2v2.𝒪est[:final] ≈ 4 atol = testTol

    # Taylor-Radau
    sim1v2 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(1), order_q = Val(2)), ref_setup)
    @test sim1v2.𝒪est[:final] ≈ 3 atol = testTol

    # Taylor-Lobatto
    sim1v3 = analyticless_test_convergence(dts, prob_stiff, ImplicitTaylor(order = Val(1), order_q = Val(3)), ref_setup)
    @test sim1v3.𝒪est[:final] ≈ 4 atol = testTol
end

include("nonlinear_convergence_tests.jl")
include("polynomial_reuse_tests.jl")
include("jacobian_reuse_tests.jl")
include("newton_recovery_tests.jl")
include("stage_accuracy_tests.jl")
include("factored_newton_tests.jl")
include("stage_formulation_tests.jl")
include("acceptance_tests.jl")
