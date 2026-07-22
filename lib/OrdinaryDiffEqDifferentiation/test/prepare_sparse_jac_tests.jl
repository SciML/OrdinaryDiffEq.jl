using OrdinaryDiffEqDifferentiation
using SparseArrays
using LinearAlgebra
using Test

using OrdinaryDiffEqDifferentiation: prepare_sparse_jac!, same_sparsity_structure

# Reference: the unconditional rebuild `prepare_sparse_jac!` replaced. Used to assert
# the guarded version is byte-identical to the old behavior on both paths.
function rebuild_reference!(J, jac_prototype)
    OrdinaryDiffEqDifferentiation.set_all_nzval!(jac_prototype, true)
    J .= true .* jac_prototype
    OrdinaryDiffEqDifferentiation.set_all_nzval!(J, false)
    return J
end

using SparseArrays: getcolptr
struct_eq(A, B) = getcolptr(A) == getcolptr(B) && rowvals(A) == rowvals(B)

@testset "prepare_sparse_jac!" begin
    # A non-full-diagonal sparse pattern (structural gap on the diagonal at row 3).
    proto = sparse([1, 2, 1, 2, 3, 1, 2], [1, 1, 2, 2, 2, 3, 3], collect(1.0:7.0), 3, 3)

    @testset "matching structure (fast path)" begin
        J = similar(proto)          # `similar` copies the structure, as build_J_W does
        @test same_sparsity_structure(J, proto)
        proto_before = copy(nonzeros(proto))

        prepare_sparse_jac!(J, proto)

        @test struct_eq(J, proto)                 # structure preserved
        @test all(iszero, nonzeros(J))            # values zeroed, ready for f.jac
        @test nonzeros(proto) == proto_before     # user's prototype NOT clobbered

        # identical to the old unconditional rebuild
        Jref = similar(proto)
        rebuild_reference!(Jref, copy(proto))
        @test struct_eq(J, Jref) && nonzeros(J) == nonzeros(Jref)
    end

    @testset "mismatched structure (rebuild path)" begin
        # J starts with a different (sparser) pattern than the prototype.
        J = sparse([1, 2], [1, 2], [9.0, 9.0], 3, 3)
        @test !same_sparsity_structure(
            J, sparse(
                [1, 2, 1, 2, 3, 1, 2],
                [1, 1, 2, 2, 2, 3, 3], collect(1.0:7.0), 3, 3
            )
        )

        p = copy(proto)
        prepare_sparse_jac!(J, p)

        @test struct_eq(J, proto)                 # rebuilt to the prototype's structure
        @test all(iszero, nonzeros(J))

        Jref = sparse([1, 2], [1, 2], [9.0, 9.0], 3, 3)
        rebuild_reference!(Jref, copy(proto))
        @test struct_eq(J, Jref) && nonzeros(J) == nonzeros(Jref)
    end
end
