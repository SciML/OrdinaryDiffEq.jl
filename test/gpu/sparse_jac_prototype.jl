using OrdinaryDiffEq
using LinearAlgebra, SparseArrays, Test
using CUDA, CUDA.CUSPARSE
using CUDSS
using OrdinaryDiffEqRosenbrock: Rosenbrock23
using OrdinaryDiffEqBDF: FBDF
using LinearSolve: LUFactorization, KrylovJL_GMRES
using SciMLBase: successful_retcode

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
