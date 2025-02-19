"""
    JVPCache{T} <: AbstractJacobianOperator{T} <: AbstractSciMLOperator{T}

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
@concrete struct JVPCache{T <: Real} <: SciMLOperators.AbstractSciMLOperator{T}
    jvp_op
    f
    du
    u
    p
    t
end

SciMLBase.isinplace(::JVPCache) = true
ArrayInterface.can_setindex(::JVPCache) = false
function ArrayInterface.restructure(y::JVPCache, x::JVPCache)
    @assert size(y)==size(x) "cannot restructure operators. ensure their sizes match."
    return x
end

function ConstructionBase.constructorof(::Type{<:JVPCache{T}}) where {T}
    return JVPCache{T}
end

Base.size(J::JVPCache) = (length(J.u), length(J.u))

function JVPCache(f::DiffEqBase.AbstractDiffEqFunction, du, u, p, t; autodiff)
    jvp_op = prepare_jvp(f, du, u, p, t, autodiff)
    return JVPCache{eltype(du)}(jvp_op, f, du, u, p, t)
end

function (op::JVPCache)(v, u, p, t)
    op.jvp_op(op.du, v, u, p, t)
    return res
end

function (op::JVPCache)(Jv, v, u, p, t)
    op.jvp_op(Jv, v, u, p, t)
    return Jv
end

Base.:*(J::JVPCache, v::AbstractArray) = J.jac_op(v, J.u, J.p, J.t)
Base.:*(J::JVPCache, v::Number) = J.jac_op(v, J.u, J.p, J.t)

function LinearAlgebra.mul!(
        Jv::AbstractArray, J::JVPCache, v::AbstractArray)
    J.jac_op(Jv, v, J.u, J.p, J.t)
    return Jv
end

# helper functions

function prepare_jvp(f::DiffEqBase.AbstractDiffEqFunction, du, u, p, t, autodiff)
    SciMLBase.has_jvp(f) && return f.jvp

    autodiff = autodiff isa AutoSparse ?  ADTypes.dense_ad(autodiff) : autodiff
    @assert DI.check_inplace(autodiff) "AD backend $(autodiff) doesn't support in-place problems."
    di_prep = DI.prepare_pushforward(
        f, du, autodiff, u, (u,), DI.Constant(p), DI.Constant(t))
    return (Jv, v, u, p, t) -> DI.pushforward!(f, du, (Jv,), di_prep,
            autodiff, u, (v, ), DI.Constant(p), DI.Constant(t))
end

function SciMLOperators.update_coefficients!(J::JVPCache, u, p, t) 
    J.u = u
    J.p = p
    J.t = t
end

