module OrdinaryDiffEqDefaultFastAlmostBandedMatricesExt

using FastAlmostBandedMatrices: AlmostBandedMatrix, bandpart, fillpart
import OrdinaryDiffEqDefault: has_structured_jac_prototype
import OrdinaryDiffEqDifferentiation: fill_stored!

has_structured_jac_prototype(A::AlmostBandedMatrix, f) =
    fillpart(A) isa StridedMatrix && has_structured_jac_prototype(bandpart(A), f)

function fill_stored!(A::AlmostBandedMatrix, v)
    fill_stored!(bandpart(A), v)
    fill!(fillpart(A), v)
    return A
end

end
