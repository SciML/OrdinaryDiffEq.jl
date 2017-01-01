immutable DiffCache{T, S}
    du::Vector{T}
    dual_du::Vector{S}
end

Base.@pure function DiffCache{chunk_size}(T, length, ::Type{Val{chunk_size}})
    DiffCache(zeros(T, length), zeros(Dual{chunk_size, T}, length))
end

Base.@pure DiffCache(u::AbstractArray) = DiffCache(eltype(u),length(u),Val{ForwardDiff.pickchunksize(length(u))})
Base.@pure DiffCache(u::AbstractArray,alg) = DiffCache(eltype(u),length(u),Val{get_chunksize(alg)})

get_du{T<:Dual}(dc::DiffCache, ::Type{T}) = dc.dual_du
get_du(dc::DiffCache, T) = dc.du

type MutableReference{T}
val::T
end

Base.@pure function autodiff_setup(f!, initial_x::Vector,cache)

    permf! = (fx, x) -> f!(x, fx)

    fx2 = copy(initial_x)
    jac_cfg = ForwardDiff.JacobianConfig{get_chunksize(cache)}(initial_x, initial_x)
    g! = (x, gx) -> ForwardDiff.jacobian!(gx, permf!, fx2, x, jac_cfg)

    fg! = (x, fx, gx) -> begin
        jac_res = DiffBase.DiffResult(fx, gx)
        ForwardDiff.jacobian!(jac_res, permf!, fx2, x, jac_cfg)
        DiffBase.value(jac_res)
    end

    return DifferentiableMultivariateFunction(f!, g!, fg!)
end

#Base.getindex(m::MutableReference)=m.val
#Base.setindex!(m::MutableReference{T}, v::{T},I...) = m.val=v

macro muladd(ex)
  esc(to_muladd(ex))
end

function to_muladd(ex)
  is_add_operation(ex) || return ex
  operands = collect(zip(
    to_muladd.((x->x.args[2]).(ex.args[2:end])),
    to_muladd.((x->x.args[3]).(ex.args[2:end]))))

  last_operation = :($(operands[end][1]) * $(operands[end][2]))

  foldr(last_operation, operands[1:end-1]) do xs, r
    :($(Base.muladd)($(xs[1]), $(xs[2]), $r))
  end
end
is_add_operation(ex::Expr) = ex.head == :call && !isempty(ex.args) && ex.args[1] == :+
is_add_operation(ex) = false
