module OrdinaryDiffEqDefaultBlockBandedMatricesExt

using BlockBandedMatrices: BandedBlockBandedMatrix, BlockBandedMatrix, blockbandwidths,
    subblockbandwidths
using BlockBandedMatrices.BlockArrays: blockisequal
using LinearAlgebra: UniformScaling, qr
using LinearSolve: GenericFactorization
using SciMLBase: isinplace
import OrdinaryDiffEqDefault: has_structured_jac_prototype, storage_linsolve
import OrdinaryDiffEqDifferentiation: fill_stored!

const BlockBandedJacobian = Union{BlockBandedMatrix, BandedBlockBandedMatrix}

fill_stored!(A::BlockBandedJacobian, v) = (fill_stored!(A.data, v); A)

function has_structured_jac_prototype(A::BlockBandedJacobian, f)
    return isinplace(f) && f.mass_matrix isa UniformScaling &&
        blockisequal(axes(A, 1), axes(A, 2)) && all(>=(0), blockbandwidths(A)) &&
        (!(A isa BandedBlockBandedMatrix) || all(>=(0), subblockbandwidths(A)))
end

function storage_linsolve(A::BlockBandedJacobian, f)
    return has_structured_jac_prototype(A, f) ? GenericFactorization(; fact_alg = qr) : nothing
end

end
