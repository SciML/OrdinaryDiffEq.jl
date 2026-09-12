using OrdinaryDiffEqDifferentiation
using SparseArrays
using SciMLOperators
using SciMLBase
using ADTypes
using Test

f!(du, u, p, t) = (du .= u; nothing)
ad_alg = AutoForwardDiff()

# `findall`/`getindex` on the mass matrix need a concrete matrix; a SciMLOperator
# (e.g. `MatrixOperator`) supports neither directly, so it must be unwrapped first.
@testset "MatrixOperator mass matrix" begin
    M = sparse([1.0 2.0; 0.0 1.0])
    mass_matrix = MatrixOperator(M; update_func = (A, u, p, t) -> M)
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; mass_matrix, jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == Matrix(M)
    @test Matrix(prob.f.jac_prototype) == Matrix(M)
end

@testset "UniformScaling mass matrix (unchanged path)" begin
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == Matrix(jac_prototype)
    @test Matrix(prob.f.jac_prototype) == Matrix(jac_prototype)
end

@testset "Plain sparse mass matrix (unchanged path)" begin
    M = sparse([1.0 2.0; 0.0 1.0])
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; mass_matrix = M, jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == Matrix(M)
    @test Matrix(prob.f.jac_prototype) == Matrix(M)
end

@testset "DiagonalOperator mass matrix" begin
    mass_matrix = DiagonalOperator([1.0, 2.0])
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; mass_matrix, jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == [1.0 0.0; 0.0 2.0]
    @test Matrix(prob.f.jac_prototype) == [1.0 0.0; 0.0 2.0]
end

# `sparsity` defaults to (and aliases) `jac_prototype`, which is how #2929's report
# constructs the problem.
@testset "MatrixOperator mass matrix, sparsity defaulted" begin
    M = sparse([1.0 2.0; 0.0 1.0])
    mass_matrix = MatrixOperator(M; update_func = (A, u, p, t) -> M)
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    odef = ODEFunction(f!; mass_matrix, jac_prototype)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.jac_prototype) == Matrix(M)
end

# Matrix-free `jac_prototype` is copied into `sparsity` by `ODEFunction`; it is not a
# sparsity pattern and must not enter the `KnownJacobianSparsityDetector` path (#4302).
@testset "Matrix-free FunctionOperator jac_prototype skips sparsity prep" begin
    using LinearAlgebra: mul!
    A = [1.0 2.0; 0.0 1.0]
    jv(v, u, p, t) = A * v
    jv(w, v, u, p, t) = mul!(w, A, v)
    Jop = FunctionOperator(jv, zeros(2), zeros(2); islinear = true)
    odef = ODEFunction(f!; jac_prototype = Jop)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    @test prob.f.sparsity === Jop
    @test OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob) === ad_alg
end

# A `DAEFunction` is fully implicit and has no `mass_matrix` field at all; the seeding must
# treat it like an identity mass matrix rather than reaching for the missing field
# (SciML/OrdinaryDiffEq.jl#1966).
@testset "DAEFunction has no mass matrix" begin
    dae_f!(res, du, u, p, t) = (res .= du .- u; nothing)
    jac_prototype = sparse([0.0 1.0; 1.0 0.0])
    daef = DAEFunction(dae_f!; jac_prototype)
    prob = DAEProblem(daef, ones(2), ones(2), (0.0, 1.0))

    prepped = OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test prepped isa ADTypes.AutoSparse
    # `DAEFunction` aliases `sparsity` to `jac_prototype`, so both gain the diagonal that
    # the `du` coefficients of the residual contribute.
    @test Matrix(prob.f.jac_prototype) == [1.0 1.0; 1.0 1.0]
end
