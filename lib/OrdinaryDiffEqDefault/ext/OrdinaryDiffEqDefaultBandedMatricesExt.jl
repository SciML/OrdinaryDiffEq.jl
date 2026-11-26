module OrdinaryDiffEqDefaultBandedMatricesExt

using BandedMatrices: BandedMatrix, bandeddata, bandwidths
using LinearAlgebra: qr
using LinearSolve: GenericFactorization
import OrdinaryDiffEqDefault: has_structured_jac_prototype, storage_linsolve
import OrdinaryDiffEqDifferentiation: fill_stored!

has_structured_jac_prototype(A::BandedMatrix{<:Any, <:StridedMatrix}, f) =
    all(>=(0), bandwidths(A))

function storage_linsolve(A::BandedMatrix, f)
    return has_structured_jac_prototype(A, f) ? GenericFactorization(; fact_alg = qr) : nothing
end

fill_stored!(A::BandedMatrix, v) = (fill!(bandeddata(A), v); A)

end
