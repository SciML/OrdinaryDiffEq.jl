"""
    JVPCache{T} <: AbstractSciMLOperator{T}

JVPCache provides a JVP operator wrapper for performing the DifferentiationInterface pushforward operation.

### Constructor

```julia
JVPCache(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t; autodiff)
```

JVPCache construction builds a DifferentiationInterface "prep" object using `prepare_pushforward!`. The "prep" object is used
when applying the operator.

### Computing the JVP

Computing the JVP is done with the DifferentiationInterface function `pushforward!`, which takes advantage of the preparation done upon construction.
"""
@concrete mutable struct JVPCache{T} <: SciMLOperators.AbstractSciMLOperator{T}
    jvp_op::Any
    f::Any
    du::Any
    u::Any
    p::Any
    t::Any
end

SciMLBase.isinplace(::JVPCache) = true
ArrayInterface.can_setindex(::JVPCache) = false
function ArrayInterface.restructure(y::JVPCache, x::JVPCache)
    @assert size(y) == size(x) "cannot restructure operators. ensure their sizes match."
    return x
end

function ConstructionBase.constructorof(::Type{<:JVPCache{T}}) where {T}
    return JVPCache{T}
end

Base.size(J::JVPCache) = (length(J.u), length(J.u))

function JVPCache(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t; autodiff)
    jvp_op = prepare_jvp(f, du, u, p, t, autodiff)
    return JVPCache{eltype(du)}(jvp_op, f, du, u, p, t)
end

function (op::JVPCache)(Jv, v, u, p, t)
    op.jvp_op(Jv, v, u, p, t)
    return Jv
end

function LinearAlgebra.mul!(
        Jv::AbstractArray, J::JVPCache, v::AbstractArray
    )
    J.jvp_op(Jv, v, J.u, J.p, J.t)
    return Jv
end

# helper functions

function prepare_jvp(f::SciMLBase.AbstractDiffEqFunction, du, u, p, t, autodiff)
    SciMLBase.has_jvp(f) && return f.jvp
    autodiff = autodiff isa AutoSparse ? ADTypes.dense_ad(autodiff) : autodiff
    @assert DI.check_inplace(autodiff) "AD backend $(autodiff) doesn't support in-place problems."
    di_prep = DI.prepare_pushforward(
        f, du, autodiff, u, (u,), DI.ConstantOrCache(p), DI.Constant(t)
    )
    return (
        Jv,
        v,
        u,
        p,
        t,
    ) -> DI.pushforward!(
        f, du, (reshape(Jv, size(du)),), di_prep, autodiff, u,
        (reshape(v, size(u)),), DI.ConstantOrCache(p), DI.Constant(t)
    )
end

function SciMLOperators.update_coefficients!(J::JVPCache, u, p, t)
    J.u = u
    J.p = p
    return J.t = t
end

function resize_JVPCache!(J::JVPCache, f, du, u, p, t, autodiff)
    J.jvp_op = prepare_jvp(f, du, u, p, t, autodiff)
    J.du = du
    return update_coefficients!(J, u, p, t)
end

function resize_JVPCache!(J::JVPCache, f, du, u, autodiff)
    J.jvp_op = prepare_jvp(f, du, u, J.p, J.t, autodiff)
    J.du = du
    return update_coefficients!(J, u, J.p, J.t)
end
