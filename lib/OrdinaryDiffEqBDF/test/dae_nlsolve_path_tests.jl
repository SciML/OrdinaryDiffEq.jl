using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, Test
using SciMLBase: DAEFunction, DAEProblem, NoInit
using ADTypes: AutoFiniteDiff
using LinearAlgebra: I
using SparseArrays

function dae_lin!(res, du, u, p, t)
    @. res = du - u
    return nothing
end
dae_lin(du, u, p, t) = @. du - u

const prob_iip = DAEProblem(dae_lin!, [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))
const prob_oop = DAEProblem{false}(dae_lin, [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))
const dae_bdf_algs = (DImplicitEuler(), DABDF2(), DFBDF(), DNordsieckBDF())

@testset "DAE BDF uses COEFFICIENT_MULTISTEP (issue #4113)" begin
    for (prob, label) in ((prob_iip, "iip"), (prob_oop, "oop"))
        @testset "$label" begin
            for alg in dae_bdf_algs
                integrator = init(
                    prob, alg; dt = 0.1, adaptive = false, save_everystep = false
                )
                step!(integrator)
                step!(integrator)
                nlsolver = integrator.cache.nlsolver
                @test nlsolver.method === OrdinaryDiffEqCore.COEFFICIENT_MULTISTEP
                @test nlsolver.cache.W_γdt ≈
                    nlsolver.α / (nlsolver.γ * integrator.dt)
            end
        end
    end
end

# A sparse matrix type without `setindex!`, as `CuSparseMatrixCSR`: `fill!` on it throws,
# so the separated DAE Jacobian `J_du` has to be zeroed through its stored values.
struct NoSetindexSparse{T} <: SparseArrays.AbstractSparseMatrix{T, Int}
    A::SparseMatrixCSC{T, Int}
end
Base.size(M::NoSetindexSparse) = size(M.A)
Base.getindex(M::NoSetindexSparse, i::Int, j::Int) = M.A[i, j]
Base.similar(M::NoSetindexSparse, ::Type{T}) where {T} = NoSetindexSparse(similar(M.A, T))
SparseArrays.nonzeros(M::NoSetindexSparse) = nonzeros(M.A)
SparseArrays.rowvals(M::NoSetindexSparse) = rowvals(M.A)
SparseArrays.nzrange(M::NoSetindexSparse, j::Integer) = nzrange(M.A, j)
SparseArrays.nnz(M::NoSetindexSparse) = nnz(M.A)

@testset "DAE J_du for a sparse jac_prototype without setindex!" begin
    sparse_dae!(res, du, u, p, t) = (res .= du .+ u; nothing)
    sparse_dae_jac!(J, du, u, p, gamma, t) = (nonzeros(J) .= 1 + gamma; nothing)
    pattern = sparse(1.0I, 4, 4)
    prototype = NoSetindexSparse(pattern)
    for (label, f) in (
            ("user jac", DAEFunction{true}(sparse_dae!; jac = sparse_dae_jac!, jac_prototype = prototype)),
            ("AD jac", DAEFunction{true}(sparse_dae!; jac_prototype = prototype)),
        )
        @testset "$label" begin
            prob = DAEProblem(f, fill(-1.0, 4), ones(4), (0.0, 1.0))
            integrator = init(prob, DFBDF(autodiff = AutoFiniteDiff()); initializealg = NoInit())
            J_du = integrator.cache.nlsolver.cache.dae_jacobians.J_du
            @test J_du isa NoSetindexSparse{Float64}
            @test rowvals(J_du) == rowvals(pattern)
            @test [nzrange(J_du, j) for j in axes(J_du, 2)] ==
                [nzrange(pattern, j) for j in axes(pattern, 2)]
            @test all(iszero, nonzeros(J_du))
        end
    end
end
