module OrdinaryDiffEqDefaultCUDAExt

using CUDA.CUSPARSE: CuSparseMatrixCSR
using LinearSolve: LinearSolve
import OrdinaryDiffEqDefault: has_structured_jac_prototype

function has_structured_jac_prototype(::CuSparseMatrixCSR, f)
    return Base.get_extension(LinearSolve, :LinearSolveCUDAExt) !== nothing &&
        Base.get_extension(LinearSolve, :LinearSolveCUDSSExt) !== nothing
end

end
