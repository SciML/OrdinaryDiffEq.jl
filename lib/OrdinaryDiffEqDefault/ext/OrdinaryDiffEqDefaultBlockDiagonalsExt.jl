module OrdinaryDiffEqDefaultBlockDiagonalsExt

using BlockDiagonals: BlockDiagonal, blocks
using LinearSolve: LinearSolve
import OrdinaryDiffEqDefault: has_structured_jac_prototype
import OrdinaryDiffEqDifferentiation: fill_stored!

@static if pkgversion(LinearSolve) >= v"5.15"
    function has_structured_jac_prototype(A::BlockDiagonal{<:Any, <:StridedMatrix}, f)
        return all(b -> size(b, 1) == size(b, 2), blocks(A))
    end
end

function fill_stored!(A::BlockDiagonal, v)
    foreach(b -> fill_stored!(b, v), blocks(A))
    return A
end

end
