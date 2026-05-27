"""
    OpaqueVoid{P, F}

A `SciMLBase.Void`-like wrapper that *de-specializes* the parameter argument
of an in-place RHS / Jacobian / time-gradient callback. The type carries the
concrete parameter type `P` so that, at call time, the
`RespecializeParams.OpaqueParams` argument is unpacked back to a `P` value
before the user's `f` is invoked.

The point is to keep the `FunctionWrappersWrapper` signature uniform on
`OpaqueParams` regardless of `P`, so a single precompiled solver path is
shared across problems whose parameter struct types differ. The unpack is a
single `unsafe_load` and is type-stable with no allocation.

`OpaqueVoid` is constructed by [`wrap_void_opaque`](@ref) during the
`AutoSpecialize` wrapping path; users do not construct it directly.
"""
struct OpaqueVoid{P, F}
    f::F
end

OpaqueVoid(::Type{P}, f::F) where {P, F} = OpaqueVoid{P, F}(f)

# 4-arg shape: rhs!(du, u, p, t), tgrad!(dT, u, p, t), jac!(J, u, p, t).
@inline function (v::OpaqueVoid{P})(a, u, op::RespecializeParams.OpaqueParams, t) where {P}
    p = RespecializeParams.unsafe_unpack(op, P)
    v.f(a, u, p, t)
    return nothing
end

# 3-arg shape: rhs!(du, u, p) for some Nonlinear in-place residuals.
@inline function (v::OpaqueVoid{P})(a, u, op::RespecializeParams.OpaqueParams) where {P}
    p = RespecializeParams.unsafe_unpack(op, P)
    v.f(a, u, p)
    return nothing
end

"""
    should_opaque_p(p) :: Bool

True if `p` is a candidate for being routed through [`OpaqueVoid`](@ref): it
must be an `isbits` value and not a `SciMLBase.NullParameters` sentinel.
Non-`isbits` and `NullParameters` payloads fall back to the existing `Void`
wrapping path unchanged.
"""
@inline should_opaque_p(p) = isbits(p) && !(p isa SciMLBase.NullParameters)

"""
    pack_p_for_opaque(p)

Return `RespecializeParams.pack(p)` if [`should_opaque_p`](@ref) returns true,
else return `p` unchanged. Single source of truth for the predicate so that
[`promote_f`](@ref) and `get_concrete_problem` stay in agreement on what
gets opaque-ified.
"""
@inline pack_p_for_opaque(p) = should_opaque_p(p) ? RespecializeParams.pack(p) : p

# Replace the third element of a Type-tuple with OpaqueParams.
@inline function _opaque_sig(::Type{Tuple{A, B, _P, T}}) where {A, B, _P, T}
    return Tuple{A, B, RespecializeParams.OpaqueParams, T}
end

"""
    wrap_void_opaque(ff, ::Type{P}, sigs::Tuple) -> FunctionWrappersWrapper

Build a `FunctionWrappersWrapper` whose signature(s) have `OpaqueParams` in
place of `typeof(p)`. `P` is the original concrete parameter type; the
produced wrapper, when invoked with the de-specialized signature, routes
through [`OpaqueVoid`](@ref) and unpacks the opaque parameter to a value of
type `P` before calling `ff`.

`sigs` is a tuple of `DataType` signatures, each of which must be a `Tuple`
type whose third element is the `p` slot.
"""
function wrap_void_opaque(ff, ::Type{P}, sigs::Tuple) where {P}
    opaque_sigs = map(_opaque_sig, sigs)
    nothings = map(_ -> Nothing, sigs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        OpaqueVoid(P, ff), opaque_sigs, nothings,
    )
end

"""
    wrapfun_iip_opaque(ff, ::Type{P}, inputs) -> FunctionWrappersWrapper

Opaque-aware analogue of [`wrapfun_iip`](@ref). `inputs` is the original
`(u0, u0, p, t)` tuple; only the third type is replaced.
"""
function wrapfun_iip_opaque(ff, ::Type{P}, inputs) where {P}
    sig = Tuple{typeof(inputs[1]), typeof(inputs[2]), typeof(inputs[3]), typeof(inputs[4])}
    return wrap_void_opaque(ff, P, (sig,))
end
