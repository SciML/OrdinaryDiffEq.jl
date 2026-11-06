using OrdinaryDiffEq
using LinearAlgebra, SparseArrays, Test
using CUDA, CUDA.CUSPARSE
using CUDSS
using OrdinaryDiffEqRosenbrock: Rosenbrock23
using OrdinaryDiffEqBDF: FBDF
using LinearSolve: LUFactorization, KrylovJL_GMRES
using SciMLBase: successful_retcode
using OrdinaryDiffEqDifferentiation: jacobian2W!, _use_allocating_sparse_W_path

#=
A stiff ODE with a sparse GPU `jac_prototype`. `jacobian2W!` forms W = J - M/dtgamma;
with the default mass matrix `I` it used to write the diagonal of `W` directly, which a
CuSparseMatrix supports neither by scalar indexing nor by broadcasting into a diagonal
view, so every implicit solver threw "Scalar indexing is disallowed" on the first step.
See https://github.com/SciML/OrdinaryDiffEq.jl/issues/1566.

The RHS is pure broadcast on purpose: `mul!` with a CuSparseMatrixCSR has no
ForwardDiff.Dual kernel and would scalar-index inside the user function, which would
fail for reasons that have nothing to do with the solver.
=#

CUDA.allowscalar(false)

const NN = 128
const DIFF = 50.0
const REACT = 100.0

rhs! = function (du, u, p, t)
    du .= DIFF .* (circshift(u, 1) .+ circshift(u, -1) .- 2 .* u) .+
        REACT .* u .* (1 .- u)
    return nothing
end

pattern = spdiagm(-1 => ones(NN - 1), 0 => ones(NN), 1 => ones(NN - 1))
pattern[1, NN] = 1.0
pattern[NN, 1] = 1.0

u0 = @. 0.5 + 0.4 * sin(2pi * (1:NN) / NN)
tspan = (0.0, 0.02)

reference = solve(
    ODEProblem(rhs!, u0, tspan), Rosenbrock23();
    save_everystep = false, abstol = 1.0e-12, reltol = 1.0e-12
).u[end]

@testset "sparse GPU jac_prototype (#1566)" begin
    # A CPU control at the same settings: if this disagrees with the reference the
    # problem setup is wrong and the GPU comparisons below mean nothing.
    cpu = solve(
        ODEProblem(ODEFunction(rhs!; jac_prototype = pattern), u0, tspan),
        Rosenbrock23(); save_everystep = false
    ).u[end]
    @test norm(cpu - reference) / norm(reference) < 1.0e-2

    if !CUDA.functional()
        @info "CUDA is not functional; skipping the GPU solves"
    else
        u0_gpu = CuArray(u0)
        proto = CuSparseMatrixCSR(CuSparseMatrixCSC(pattern))

        # The W build is what regressed, so cover the direct and the Krylov paths as
        # well as a second solver family.
        algs = (
            ("default", Rosenbrock23()),
            ("LU", Rosenbrock23(linsolve = LUFactorization())),
            ("GMRES", Rosenbrock23(linsolve = KrylovJL_GMRES(), concrete_jac = true)),
            ("FBDF", FBDF(linsolve = KrylovJL_GMRES(), concrete_jac = true)),
        )

        @testset "$name" for (name, alg) in algs
            sol = solve(
                ODEProblem(ODEFunction(rhs!; jac_prototype = proto), u0_gpu, tspan),
                alg; save_everystep = false
            )
            @test successful_retcode(sol)
            @test norm(Array(sol.u[end]) - reference) / norm(reference) < 1.0e-2
        end
    end
end

# A prototype whose diagonal is not fully stored: `W = J - M/dtgamma` fills those
# entries in, so the sum has more stored entries than the preallocated `W`. Both
# patterns used above store the whole diagonal, so this case needs its own test.
@testset "prototype with a structural zero on the diagonal (#4452)" begin
    m = 6
    Jcpu = sparse(
        [1, 2, 3, 4, 5, 6, 1, 2], [1, 2, 4, 4, 5, 6, 2, 1],
        [2.0, 3.0, 1.0, 4.0, 5.0, 6.0, 0.5, 0.7], m, m
    )
    stored = count(i -> Jcpu[i, i] != 0, 1:m)
    @test stored < m          # the point of the test is that some are missing

    dtgamma = 0.25
    reference = Matrix(Jcpu) - inv(dtgamma) * Matrix(I, m, m)

    Jgpu = CUSPARSE.CuSparseMatrixCSR(CUSPARSE.CuSparseMatrixCSC(Jcpu))
    W = CUSPARSE.CuSparseMatrixCSR(CUSPARSE.CuSparseMatrixCSC(copy(Jcpu)))
    @test _use_allocating_sparse_W_path(W)

    jacobian2W!(W, I, dtgamma, Jgpu)

    # Reading back is allowed to index scalars; it is the check, not the code
    # under test. What matters is that the filled-in diagonal is present and right.
    got = CUDA.@allowscalar Matrix(CUSPARSE.CuSparseMatrixCSC(W))
    @test nnz(W) == nnz(Jcpu) + (m - stored)
    @test got ≈ reference
end
