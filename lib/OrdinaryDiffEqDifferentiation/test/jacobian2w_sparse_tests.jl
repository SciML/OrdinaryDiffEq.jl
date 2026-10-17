using OrdinaryDiffEqDifferentiation, OrdinaryDiffEqCore
using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqBDF
using LinearAlgebra
using SparseArrays
using Test

# `jacobian2W!` forms W = J - M/dtgamma. With a scalar mass matrix it writes the
# diagonal, which a sparse GPU `W` supports neither by scalar indexing nor by
# broadcasting into a diagonal view, so that case takes the allocating path instead.
# The allocating path is only selected for GPU storage, so the equivalence of the two
# formulas is what can be checked on CPU; the GPU dispatch itself is covered by
# test/gpu/sparse_jac_prototype.jl.

J = sparse([3.0 1.0; 0.5 2.0])
dtgamma = 0.25
invdtgamma = inv(dtgamma)

for mm in (I, 2.0 * I)
    λ = -OrdinaryDiffEqDifferentiation._scalar_massmatrix_λ(mm)

    W = similar(J)
    fill!(nonzeros(W), 0)
    OrdinaryDiffEqDifferentiation.jacobian2W!(W, mm, dtgamma, J)

    # The formula the sparse-GPU branch uses must agree with the diagonal write.
    allocating = J + (λ * invdtgamma) * I
    @test Matrix(W) ≈ Matrix(allocating)
    @test Matrix(W) ≈ Matrix(J) - OrdinaryDiffEqDifferentiation._scalar_massmatrix_λ(mm) *
        invdtgamma * Matrix(I, 2, 2)
end

# Bitwise check: the pattern-stable fill must produce exactly the values the
# broadcast / copyto! paths produced, including a stored -0.0 in J.
function bits(A)
    return map(bitstring, Matrix(A))
end

Jzero = sparse([1, 1, 2, 2], [1, 2, 1, 2], [3.0, -0.0, 0.5, 2.0])
for mm in (I, 3.5 * I, Diagonal([1.0, 0.0]))
    W = similar(Jzero)
    OrdinaryDiffEqDifferentiation.jacobian2W!(W, mm, dtgamma, Jzero)
    # Reference semantics: scalar `M` copies `J` bitwise and shifts the
    # diagonal; `Diagonal` `M` broadcasts `muladd(-M, invdtgamma, J)` with
    # `(-M)[i,j] = +0.0` off-diagonal.
    Wref = if mm isa UniformScaling
        R = Matrix(Jzero)
        λ = -mm.λ
        for i in 1:2
            R[i, i] = muladd(λ, invdtgamma, Jzero[i, i])
        end
        R
    else
        muladd.(Matrix(-mm), invdtgamma, Matrix(Jzero))
    end
    @test bits(W) == bits(Wref)
end

# Widening: W is grown once to `W ∪ J ∪ nz(M)`, then the pattern stays fixed.
Jmissing = sparse([1, 2], [2, 1], [2.0, 3.0], 2, 2)
Wmissing = copy(Jmissing)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wmissing, I, dtgamma, Jmissing)
@test Matrix(Wmissing) ≈ Matrix(Jmissing) - invdtgamma * Matrix(I, 2, 2)
@test nnz(Wmissing) == nnz(Jmissing) + 2

# Partial diagonal: the single missing diagonal slot is inserted.
Jpartial = sparse([1, 1, 2], [1, 2, 1], [3.0, 1.0, 0.5], 2, 2)
Wpartial = copy(Jpartial)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wpartial, I, dtgamma, Jpartial)
@test Matrix(Wpartial) ≈ Matrix(Jpartial) - invdtgamma * Matrix(I, 2, 2)
@test nnz(Wpartial) == nnz(Jpartial) + 1

# A `Diagonal` mass matrix only needs its nonzero diagonal slots (nz(M)).
Jdiag = sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0, 2.0, 0.5, 0.5], 2, 2)
Wdiag = copy(Jdiag)
OrdinaryDiffEqDifferentiation.jacobian2W!(
    Wdiag, Diagonal([2.0, 0.0]), dtgamma, Jdiag
)
@test Matrix(Wdiag) ≈ Matrix(Jdiag) - invdtgamma * Diagonal([2.0, 0.0])
@test nnz(Wdiag) == nnz(Jdiag)

# A mass matrix that is neither scalar nor Diagonal keeps the broadcast path.
Mfull = sparse([1.0 0.1; 0.0 1.0])
Wfull = similar(J)
fill!(nonzeros(Wfull), 0)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wfull, Mfull, dtgamma, J)
@test Matrix(Wfull) ≈ Matrix(J) - invdtgamma * Matrix(Mfull)

# Pattern stability: with a Jacobian entry that is exactly 0.0 in some calls,
# the broadcast path would drop the slot and `nnz(W)` would oscillate. The
# pattern-stable path keeps every stored slot.
Jvarying = sparse([1, 1, 2, 2], [1, 2, 1, 2], [1.0, 1.0, 0.5, -1.0])
Wstable = similar(Jvarying)
nnz_hist = Int[]
for val in (0.0, 1.0, 0.0, -2.0)
    Jvarying[1, 2] = val
    OrdinaryDiffEqDifferentiation.jacobian2W!(
        Wstable, Diagonal([1.0, 0.0]), dtgamma, Jvarying
    )
    push!(nnz_hist, nnz(Wstable))
end
@test allequal(nnz_hist)
@test Matrix(Wstable) ≈ Matrix(Jvarying) - invdtgamma * Diagonal([1.0, 0.0])

# CPU sparse storage is scalar-indexable, so it keeps the in-place path; a
# dense matrix is not sparse at all. Only GPU storage takes the allocating branch.
W = similar(J)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(W)
@test !OrdinaryDiffEqDifferentiation._use_allocating_sparse_W_path(ones(2, 2))

# The dense `Matrix` method must be unaffected by the sparse guard.
Jd = Matrix(J)
Wd = similar(Jd)
OrdinaryDiffEqDifferentiation.jacobian2W!(Wd, I, dtgamma, Jd)
@test Wd ≈ Jd - invdtgamma * Matrix(I, 2, 2)

@testset "nnz(W) constant across steps of a mass-matrix DAE" begin
    # Mass-matrix DAE: du1 = u2, du2 = -u1 + c(t) u3, 0 = u1 + u2 - u3.
    # The (2,3) Jacobian entry c(t) is exactly 0.0 for t < 0.5, so the
    # broadcast path drops it from W's pattern for part of the solve.
    function dae_rhs!(du, u, p, t)
        du[1] = u[2]
        du[2] = -u[1] + (t < 0.5 ? zero(t) : one(t)) * u[3]
        du[3] = u[1] + u[2] - u[3]
        return nothing
    end
    function dae_jac!(J, u, p, t)
        J[1, 2] = 1.0
        J[2, 1] = -1.0
        J[2, 3] = t < 0.5 ? 0.0 : 1.0
        J[3, 1] = 1.0
        J[3, 2] = 1.0
        J[3, 3] = -1.0
        return nothing
    end
    jac_prototype = sparse(
        [1, 2, 2, 3, 3, 3], [2, 1, 3, 1, 2, 3], ones(6), 3, 3
    )
    f = ODEFunction(
        dae_rhs!; jac = dae_jac!, jac_prototype = jac_prototype,
        mass_matrix = Diagonal([1.0, 1.0, 0.0])
    )
    prob = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0))

    function nnz_over_steps(alg, get_W)
        integrator = init(prob, alg; abstol = 1.0e-8, reltol = 1.0e-8)
        nnzs = Int[]
        for _ in integrator
            push!(nnzs, nnz(get_W(integrator)))
        end
        return nnzs
    end

    # `always_new` forces FBDF to rebuild W on every nonlinear solve so its
    # sparse broadcast (or lack of one) is actually exercised.
    for (name, alg, get_W) in (
            ("Rodas5P", Rodas5P(), integ -> integ.cache.W),
            (
                "FBDF",
                FBDF(;
                    nlsolve = OrdinaryDiffEqBDF.NLNewton(always_new = true)
                ),
                integ -> integ.cache.nlsolver.cache.W,
            ),
        )
        nnzs = nnz_over_steps(alg, get_W)
        @testset "$name" begin
            @test length(nnzs) > 1
            @test allequal(nnzs)
        end
    end
end
