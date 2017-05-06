immutable DiffEqDiffTools.DiffCache{T, S}
    du::Vector{T}
    dual_du::Vector{S}
end

Base.@pure function DiffEqDiffTools.DiffCache{chunk_size}(T, length, ::Type{Val{chunk_size}})
    DiffEqDiffTools.DiffCache(zeros(T, length), zeros(Dual{chunk_size, T}, length))
end

Base.@pure DiffEqDiffTools.DiffCache(u::AbstractArray) = DiffEqDiffTools.DiffCache(eltype(u),length(u),Val{ForwardDiff.pickchunksize(length(u))})
Base.@pure DiffEqDiffTools.DiffCache(u::AbstractArray,nlsolve) = DiffEqDiffTools.DiffCache(eltype(u),length(u),Val{DiffEqDiffTools.get_chunksize(nlsolve)})
Base.@pure DiffEqDiffTools.DiffCache{CS}(u::AbstractArray,T::Type{Val{CS}}) = DiffEqDiffTools.DiffCache(eltype(u),length(u),T)

DiffEqDiffTools.get_du{T<:Dual}(dc::DiffEqDiffTools.DiffCache, ::Type{T}) = dc.dual_du
DiffEqDiffTools.get_du(dc::DiffEqDiffTools.DiffCache, T) = dc.du

DiffEqDiffTools.realtype{T}(::Type{T}) = T
DiffEqDiffTools.realtype{T}(::Type{Complex{T}}) = T

# Default nlsolve behavior, should move to DiffEqDiffTools.jl

Base.@pure DiffEqDiffTools.determine_chunksize(u,alg::DEAlgorithm) = DiffEqDiffTools.determine_chunksize(u,DiffEqDiffTools.get_chunksize(alg))
Base.@pure function DiffEqDiffTools.determine_chunksize(u,CS)
  if CS != 0
    return CS
  else
    return ForwardDiff.pickchunksize(length(u))
  end
end

function DiffEqDiffTools.autodiff_setup{CS}(f!, initial_x::Vector,chunk_size::Type{Val{CS}})

    permf! = (fx, x) -> f!(x, fx)

    fx2 = copy(initial_x)
    jac_cfg = ForwardDiff.JacobianConfig{CS}(initial_x, initial_x)
    g! = (x, gx) -> ForwardDiff.jacobian!(gx, permf!, fx2, x, jac_cfg)

    fg! = (x, fx, gx) -> begin
        jac_res = DiffBase.DiffResult(fx, gx)
        ForwardDiff.jacobian!(jac_res, permf!, fx2, x, jac_cfg)
        DiffBase.value(jac_res)
    end

    return DifferentiableMultivariateFunction(f!, g!, fg!)
end

function non_DiffEqDiffTools.autodiff_setup(f!, initial_x::Vector)
  DifferentiableMultivariateFunction(f!)
end

immutable DiffEqDiffTools.NLSOLVEJL_SETUP{CS,AD} end
Base.@pure DiffEqDiffTools.NLSOLVEJL_SETUP(;chunk_size=0,autodiff=true) = DiffEqDiffTools.NLSOLVEJL_SETUP{chunk_size,autodiff}()
(p::DiffEqDiffTools.NLSOLVEJL_SETUP)(f,u0) = (res=NLsolve.nlsolve(f,u0); res.zero)
function (p::DiffEqDiffTools.NLSOLVEJL_SETUP{CS,AD}){CS,AD}(::Type{Val{:init}},f,u0_prototype)
  if AD
    return non_DiffEqDiffTools.autodiff_setup(f,u0_prototype)
  else
    return DiffEqDiffTools.autodiff_setup(f,u0_prototype,Val{DiffEqDiffTools.determine_chunksize(initial_x,CS)})
  end
end

DiffEqDiffTools.get_chunksize(x) = 0
DiffEqDiffTools.get_chunksize{CS,AD}(x::DiffEqDiffTools.NLSOLVEJL_SETUP{CS,AD}) = CS

export DiffEqDiffTools.NLSOLVEJL_SETUP
