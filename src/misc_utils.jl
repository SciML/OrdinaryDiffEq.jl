immutable DiffCache{T, S}
    du::Vector{T}
    dual_du::Vector{S}
end

Base.@pure function DiffCache{chunk_size}(T, length, ::Type{Val{chunk_size}})
    DiffCache(zeros(T, length), zeros(Dual{chunk_size, T}, length))
end

Base.@pure DiffCache(u::AbstractArray) = DiffCache(eltype(u),length(u),Val{ForwardDiff.pickchunksize(length(u))})
Base.@pure DiffCache(u::AbstractArray,alg) = DiffCache(eltype(u),length(u),Val{get_chunksize(alg)})
Base.@pure DiffCache{CS}(u::AbstractArray,T::Type{Val{CS}}) = DiffCache(eltype(u),length(u),T)

get_du{T<:Dual}(dc::DiffCache, ::Type{T}) = dc.dual_du
get_du(dc::DiffCache, T) = dc.du

function autodiff_setup(f!, initial_x::Vector,alg)
  autodiff_setup(f!, initial_x,Val{determine_chunksize(initial_x,alg)})
end

function non_autodiff_setup(f!, initial_x::Vector,alg)
  non_autodiff_setup(f!, initial_x,Val{determine_chunksize(initial_x,alg)})
end

Base.@pure function determine_chunksize(u,alg)
  if get_chunksize(alg) != 0
    return get_chunksize(alg)
  else
    return ForwardDiff.pickchunksize(length(u))
  end
end

function autodiff_setup{CS}(f!, initial_x::Vector,chunk_size::Type{Val{CS}})

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

function non_autodiff_setup{CS}(f!, initial_x::Vector,chunk_size::Type{Val{CS}})
  DifferentiableMultivariateFunction(f!)
end

realtype{T}(::Type{T}) = T
realtype{T}(::Type{Complex{T}}) = T
