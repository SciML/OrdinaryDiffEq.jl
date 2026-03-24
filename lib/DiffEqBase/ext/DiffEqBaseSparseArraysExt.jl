module DiffEqBaseSparseArraysExt

import SparseArrays
import DiffEqBase: NAN_CHECK, INFINITE_OR_GIANT

function NAN_CHECK(x::SparseArrays.AbstractSparseMatrixCSC)
    return any(NAN_CHECK, SparseArrays.nonzeros(x))
end
function INFINITE_OR_GIANT(x::SparseArrays.AbstractSparseMatrixCSC)
    return any(INFINITE_OR_GIANT, SparseArrays.nonzeros(x))
end

end
