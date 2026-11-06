using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore
using SciMLOperators
using SciMLBase
using LinearAlgebra
using Test

f!(du, u, p, t) = (du[1] = -(1 + t) * u[1]; du[2] = -2u[2]; nothing)
update_jac!(M, u, p, t) = (M[1, 1] = -(1 + t); M[2, 2] = -2.0; nothing)

operator_problem() = ODEProblem(
    ODEFunction(
        f!; jac_prototype = MatrixOperator([-1.0 0.0; 0.0 -2.0]; update_func! = update_jac!)
    ),
    [1.0, 1.0], (0.0, 1.0)
)

@testset "an operator Jacobian reports itself stale after it moves" begin
    integrator = init(operator_problem(), ImplicitEuler(); dt = 0.1, adaptive = false)
    step!(integrator)
    W = OrdinaryDiffEqCore.get_W(integrator.cache.nlsolver)
    @test W isa SciMLOperators.WOperator
    @test W.J isa SciMLOperators.AbstractSciMLOperator

    SciMLOperators.mark_jacobian_current!(W)
    @test !SciMLOperators.jacobian_stale(W)

    step!(integrator)
    @test SciMLOperators.jacobian_stale(W)
end

@testset "a concrete Jacobian still integrates correctly" begin
    prob = ODEProblem(f!, [1.0, 1.0], (0.0, 1.0))
    reference = solve(prob, ImplicitEuler(); dt = 0.001, adaptive = false)
    coarse = solve(prob, ImplicitEuler(); dt = 0.05, adaptive = false)
    @test SciMLBase.successful_retcode(coarse)
    @test maximum(abs.(coarse.u[end] .- reference.u[end])) < 0.05
end

@testset "an operator Jacobian still integrates correctly" begin
    exact = [exp(-1.5), exp(-2.0)]
    operator = solve(operator_problem(), ImplicitEuler(); dt = 0.001, adaptive = false)
    dense = solve(
        ODEProblem(f!, [1.0, 1.0], (0.0, 1.0)), ImplicitEuler();
        dt = 0.001, adaptive = false
    )
    @test SciMLBase.successful_retcode(operator)
    @test maximum(abs.(operator.u[end] .- exact)) < 1.0e-3
    @test maximum(abs.(dense.u[end] .- exact)) < 1.0e-3
    @test operator.u[end] ≈ dense.u[end] rtol = 1.0e-4
end
