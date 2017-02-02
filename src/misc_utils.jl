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

type MutableReference{T}
val::T
end

Base.@pure function autodiff_setup(f!, initial_x::Vector,alg)
  autodiff_setup(f!, initial_x,Val{determine_chunksize(initial_x,alg)})
end

Base.@pure function determine_chunksize(u,alg)
  if get_chunksize(alg) != 0
    return get_chunksize(alg)
  else
    return ForwardDiff.pickchunksize(length(u))
  end
end

Base.@pure function autodiff_setup{CS}(f!, initial_x::Vector,chunk_size::Type{Val{CS}})

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

#Base.getindex(m::MutableReference)=m.val
#Base.setindex!(m::MutableReference{T}, v::{T},I...) = m.val=v

macro muladd(ex)
  esc(to_muladd(ex))
end

function to_muladd(ex)
  is_add_operation(ex) || return ex

  all_operands = ex.args[2:end]
  mul_operands = filter(is_mul_operation, all_operands)
  odd_operands = filter(x->!is_mul_operation(x), all_operands)

  muladd_operands = collect(zip(
    to_muladd.((x->x.args[2]).(mul_operands)),
    to_muladd.((x->x.args[3]).(mul_operands))))

  if isempty(odd_operands)
    to_be_muladded = muladd_operands[1:end-1]
    last_operation = :($(muladd_operands[end][1]) * $(muladd_operands[end][2]))
  else
    to_be_muladded = muladd_operands
    last_operation = make_addition(odd_operands)
  end

  foldr(last_operation, to_be_muladded) do xs, r
    :($(Base.muladd)($(xs[1]), $(xs[2]), $r))
  end
end

is_operation(ex::Expr, op::Symbol) = ex.head == :call && !isempty(ex.args) && ex.args[1] == op
is_operation(ex, op::Symbol) = false

is_add_operation(ex) = is_operation(ex, :+)
is_mul_operation(ex) = is_operation(ex, :*)

make_addition(args) = length(args) == 1 ? args[1] : Expr(:call, :+, args...)

realtype{T}(::Type{T}) = T
realtype{T}(::Type{Complex{T}}) = T
