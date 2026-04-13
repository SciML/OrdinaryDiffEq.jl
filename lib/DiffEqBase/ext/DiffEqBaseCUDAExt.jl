module DiffEqBaseCUDAExt

using DiffEqBase, CUDA

function DiffEqBase.ODE_DEFAULT_NORM(
        u::CuArray{T}, t
    ) where {T <: Union{AbstractFloat, Complex}}
    return sqrt(sum(DiffEqBase.sse, u; init = DiffEqBase.sse(zero(T))) / DiffEqBase.totallength(u))
end

end
