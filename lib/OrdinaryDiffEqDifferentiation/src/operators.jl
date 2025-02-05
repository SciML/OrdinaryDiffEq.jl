abstract type AbstractJacobianOperator{T} <: SciMLOperators.AbstractSciMLOperator{T} end

ArrayInterface.can_setindex(::AbstractJacobianOperator) = false
function ArrayInterface.restructure(
        y::AbstractJacobianOperator, x::AbstractJacobianOperator
)
    @assert size(y)==size(x) "cannot restructure operators. ensure their sizes match."
    return x
end

abstract type AbstractMode end

struct VJP <: AbstractMode end
struct JVP <: AbstractMode end

flip_mode(::VJP) = JVP()
flip_mode(::JVP) = VJP()

"""
    JacobianOperator{iip, T} <: AbstractJacobianOperator{T} <: AbstractSciMLOperator{T}

A Jacobian Operator Provides both JVP and VJP without materializing either (if possible).

### Constructor

```julia
JacobianOperator(prob::AbstractNonlinearProblem, fu, u; jvp_autodiff = nothing,
    vjp_autodiff = nothing, skip_vjp::Val = Val(false), skip_jvp::Val = Val(false))
```

By default, the `JacobianOperator` will compute `JVP`. Use `Base.adjoint` or
`Base.transpose` to switch to `VJP`.

### Computing the VJP

Computing the VJP is done according to the following rules:

  - If `f` has a `vjp` method, then we use that.
  - If `f` has a `jac` method and no `vjp_autodiff` is provided, then we use `jac * v`.
  - If `vjp_autodiff` is provided we using DifferentiationInterface.jl to compute the VJP.

### Computing the JVP

Computing the JVP is done according to the following rules:

  - If `f` has a `jvp` method, then we use that.
  - If `f` has a `jac` method and no `jvp_autodiff` is provided, then we use `v * jac`.
  - If `jvp_autodiff` is provided we using DifferentiationInterface.jl to compute the JVP.

### Special Case (Number)

For Number inputs, VJP and JVP are not distinct. Hence, if either `vjp` or `jvp` is
provided, then we use that. If neither is provided, then we use `v * jac` if `jac` is
provided. Finally, we use the respective autodiff methods to compute the derivative
using DifferentiationInterface.jl and multiply by `v`.

### Methods Provided

!!! warning

    Currently it is expected that `p` during problem construction is same as `p` during
    operator evaluation. This restriction will be lifted in the future.

  - `(op::JacobianOperator)(v, u, p)`: Computes `∂f(u, p)/∂u * v` or `∂f(u, p)/∂uᵀ * v`.
  - `(op::JacobianOperator)(res, v, u, p)`: Computes `∂f(u, p)/∂u * v` or `∂f(u, p)/∂uᵀ * v`
    and stores the result in `res`.

See also [`VecJacOperator`](@ref) and [`JacVecOperator`](@ref).
"""
@concrete struct JacobianOperator{iip, T <: Real} <: AbstractJacobianOperator{T}
    mode <: AbstractMode

    jvp_op::Any
    vjp_op::Any

    size::Any

    output_cache::Any
    input_cache::Any
end

SciMLBase.isinplace(::JacobianOperator{iip}) where {iip} = iip

function ConstructionBase.constructorof(::Type{<:JacobianOperator{iip, T}}) where {iip, T}
    return JacobianOperator{iip, T}
end

Base.size(J::JacobianOperator) = J.size
Base.size(J::JacobianOperator, d::Integer) = J.size[d]

for op in (:adjoint, :transpose)
    @eval function Base.$(op)(operator::JacobianOperator{iip, T}) where {iip, T}
        return JacobianOperator{iip, T}(
            flip_mode(operator.mode), operator.jvp_op, operator.vjp_op,
            reverse(operator.size), operator.input_cache, operator.output_cache)
    end
end

function JacobianOperator(f::DiffEqBase.AbstractDiffEqFunction, u, p, t, fu = nothing; jvp_autodiff = nothing,
        vjp_autodiff = nothing, skip_vjp::Val = Val(false), skip_jvp::Val = Val(false))

    @assert !(skip_vjp === Val(true) && skip_jvp === Val(true)) "Cannot skip both vjp and jvp \
                                                       construction."

    isnothing(fu) ? (fu = !SciMLBase.isinplace(f) ? f(u, p, t) : u) : fu

    iip = SciMLBase.isinplace(f)
    T = promote_type(eltype(u), eltype(fu))

    vjp_autodiff = vjp_autodiff
    vjp_op = prepare_vjp(skip_vjp, f, u, p, t, fu; autodiff = vjp_autodiff)

    jvp_autodiff = jvp_autodiff
    jvp_op = prepare_jvp(skip_jvp, f, u, p, t, fu; autodiff = jvp_autodiff)

    output_cache = fu isa Number ? T(fu) : similar(fu, T)
    input_cache = u isa Number ? T(u) : similar(u, T)

    return JacobianOperator{iip, T}(
        JVP(), jvp_op, vjp_op, (length(fu), length(u)), output_cache, input_cache)
end

function (op::JacobianOperator)(v, u, p, t)
    if op.mode isa VJP
        if SciMLBase.isinplace(op)
            res = zero(op.output_cache)
            op.vjp_op(res, v, u, p, t)
            return res
        end
        return op.vjp_op(v, u, p, t)
    else
        if SciMLBase.isinplace(op)
            res = zero(op.output_cache)
            op.jvp_op(res, v, u, p, t)
            return res
        end
        return op.jvp_op(v, u, p, t)
    end
end

function (op::JacobianOperator)(::Number, ::Number, _, __)
    error("Inplace Jacobian Operator not possible for scalars.")
end

function (op::JacobianOperator)(Jv, v, u, p, t)
    if op.mode isa VJP
        if SciMLBase.isinplace(op)
            op.vjp_op(Jv, v, u, p, t)
        else
            copyto!(Jv, op.vjp_op(v, u, p, t))
        end
    else
        if SciMLBase.isinplace(op)
            op.jvp_op(Jv, v, u, p, t)
        else
            copyto!(Jv, op.jvp_op(v, u, p, t))
        end
    end
    return Jv
end

"""
    StatefulJacobianOperator(jac_op::JacobianOperator, u, p, t)

Wrapper over a [`JacobianOperator`](@ref) which stores the input `u`, `p` and `t`, and defines
`mul!` and `*` for computing VJPs and JVPs.
"""
@concrete mutable struct StatefulJacobianOperator{M <: AbstractMode, T} <:
                 AbstractJacobianOperator{T}
    mode::M
    jac_op <: JacobianOperator
    u
    p
    t

    function StatefulJacobianOperator(jac_op::JacobianOperator, u, p, t)
        return new{
            typeof(jac_op.mode), eltype(jac_op), typeof(jac_op), typeof(u), typeof(p), typeof(t)}(
            jac_op.mode, jac_op, u, p, t)
    end
end

Base.size(J::StatefulJacobianOperator) = size(J.jac_op)
Base.size(J::StatefulJacobianOperator, d::Integer) = size(J.jac_op, d)

for op in (:adjoint, :transpose)
    @eval function Base.$(op)(operator::StatefulJacobianOperator)
        return StatefulJacobianOperator($(op)(operator.jac_op), operator.u, operator.p, operator.t)
    end
end

Base.:*(J::StatefulJacobianOperator, v::AbstractArray) = J.jac_op(v, J.u, J.p, J.t)
Base.:*(J::StatefulJacobianOperator, v::Number) = J.jac_op(v, J.u, J.p, J.t)

function LinearAlgebra.mul!(
        Jv::AbstractArray, J::StatefulJacobianOperator, v::AbstractArray)
    J.jac_op(Jv, v, J.u, J.p, J.t)
    return Jv
end



# helper functions

prepare_vjp(::Val{true}, args...; kwargs...) = nothing

function prepare_vjp(
        ::Val{false}, f::DiffEqBase.AbstractDiffEqFunction, u, p, t, fu; autodiff = nothing)
    SciMLBase.has_vjp(f) && return f.vjp

    autodiff isa AutoSparse ? autodiff = ADTypes.dense_ad(autodiff) : autodiff = autodiff

    if isnothing(autodiff) && SciMLBase.has_jac(f)
        if SciMLBase.isinplace(f)
            jac_cache = similar(u, eltype(fu), length(fu), length(u))
            return @closure (vJ, v, u, p, t) -> begin
                f.jac(jac_cache, u, p, t)
                LinearAlgebra.mul!(vec(vJ), jac_cache', vec(v))
                return
            end
            return vjp_op
        else
            return @closure (v, u, p, t) -> reshape(f.jac(u, p, t)' * vec(v), size(u))
        end
    end

    @assert autodiff!==nothing "`vjp_autodiff` must be provided if `f` doesn't have \
                               analytic `vjp` or `jac`."

    if SciMLBase.isinplace(f)
        @assert DI.check_inplace(autodiff) "AD backend $(autodiff) doesn't support in-place problems."

        fu_cache = copy(fu)

        di_prep = DI.prepare_pullback(
            f, fu_cache, autodiff, u, (fu,), DI.Constant(p), DI.Constant(t))
        return @closure (vJ, v, u, p, t) -> begin
            DI.pullback!(f, fu_cache, (reshape(vJ, size(u)),), di_prep, autodiff, u,
                (reshape(v, size(fu_cache)),), DI.Constant(p), DI.Constant(t))
            return
        end
    else
        di_prep = DI.prepare_pullback(f, autodiff, u, (fu,), DI.Constant(p), DI.Constant(t))
        return @closure (v, u, p, t) -> begin
            return only(DI.pullback(
                f, di_prep, autodiff, u, (reshape(v, size(fu)),), DI.Constant(p), DI.Constant(t)))
        end
    end
end


prepare_jvp(skip::Val{true}, args...; kwargs...) = nothing

function prepare_jvp(
        ::Val{false}, f::DiffEqBase.AbstractDiffEqFunction, u, p, t, fu; autodiff = nothing)

    SciMLBase.has_jvp(f) && return f.jvp

    autodiff isa AutoSparse ? autodiff = ADTypes.dense_ad(autodiff) : autodiff = autodiff

    if isnothing(autodiff) && SciMLBase.has_jac(f)
        if SciMLBase.isinplace(f)
            jac_cache = similar(u, eltype(fu), length(fu), length(u))
            return @closure (Jv, v, u, p, t) -> begin
                f.jac(jac_cache, u, p, t)
                LinearAlgebra.mul!(vec(Jv), jac_cache, vec(v))'
                return
            end
        else
            return @closure (v, u, p, t) -> reshape(f.jac(u, p, t) * vec(v), size(u))
        end
    end

    @assert autodiff!==nothing "`jvp_autodiff` must be provided if `f` doesn't have \
                                analytic `jvp` or `jac`"

    if SciMLBase.isinplace(f)
        @assert DI.check_inplace(autodiff) "AD backend $(autodiff) doesn't support in-place problems."

        fu_cache = copy(fu)
        di_prep = DI.prepare_pushforward(
            f, fu_cache, autodiff, u, (u,), DI.Constant(p), DI.Constant(t))
        return (Jv, v, u, p, t) -> begin
            return DI.pushforward!(f, fu_cache, (reshape(Jv, size(fu_cache)),), di_prep,
                autodiff, u, (reshape(v, size(u)),), DI.Constant(p), DI.Constant(t))
        end
    else
        di_prep = DI.prepare_pushforward(f, autodiff, u, (u,), DI.Constant(p), DI.Constant(t))
        return @closure (v, u, p, t) -> begin
            return only(DI.pushforward(
                f, di_prep, autodiff, u, (reshape(v, size(u)),), DI.Constant(p), DI.Constant(t)))
        end
    end
end

function SciMLOperators.update_coefficients!(J::StatefulJacobianOperator, u, p, t) 
    J.u = u
    J.p = p
    J.t = t
end

