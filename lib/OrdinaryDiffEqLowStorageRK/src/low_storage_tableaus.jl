@enum LowStorageRKForm TwoN TwoC

struct LowStorageRKTableau{form, N, T, T2} <: OrdinaryDiffEqConstantCache
    A2end::NTuple{N, T}
    B1::T
    B2end::NTuple{N, T}
    c2end::NTuple{N, T2}
end

function LowStorageRKTableau{form}(
        A2end::NTuple{N, T}, B1::T, B2end::NTuple{N, T}, c2end::NTuple{N, T2}
    ) where {form, N, T, T2}
    return LowStorageRKTableau{form, N, T, T2}(A2end, B1, B2end, c2end)
end
