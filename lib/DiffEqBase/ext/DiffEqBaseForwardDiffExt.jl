module DiffEqBaseForwardDiffExt

using DiffEqBase, ForwardDiff
using DiffEqBase.ArrayInterface
using DiffEqBase: Void, FunctionWrappersWrappers, OrdinaryDiffEqTag,
    AbstractTimeseriesSolution,
    RecursiveArrayTools, reduce_tup, _promote_tspan, has_continuous_callback
import DiffEqBase: hasdualpromote, wrapfun_oop, wrapfun_iip, prob2dtmin,
    promote_tspan, ODE_DEFAULT_NORM
import SciMLBase: isdualtype, DualEltypeChecker, sse, __sum

const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}
dualgen(::Type{T}) where {T} = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, T}, T, 1}
dualgen(::Type{T}, ::Val{CS}) where {T, CS} = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, T}, T, CS}

const NORECOMPILE_IIP_SUPPORTED_ARGS = (
    Tuple{
        Vector{Float64}, Vector{Float64},
        Vector{Float64}, Float64,
    },
    Tuple{
        Vector{Float64}, Vector{Float64},
        SciMLBase.NullParameters, Float64,
    },
)

const oop_arglists = (
    Tuple{Vector{Float64}, Vector{Float64}, Float64},
    Tuple{Vector{Float64}, SciMLBase.NullParameters, Float64},
    Tuple{Vector{Float64}, Vector{Float64}, dualT},
    Tuple{Vector{dualT}, Vector{Float64}, Float64},
    Tuple{Vector{dualT}, SciMLBase.NullParameters, Float64},
    Tuple{Vector{Float64}, SciMLBase.NullParameters, dualT},
)

const NORECOMPILE_OOP_SUPPORTED_ARGS = (
    Tuple{
        Vector{Float64},
        Vector{Float64}, Float64,
    },
    Tuple{
        Vector{Float64},
        SciMLBase.NullParameters, Float64,
    },
)
const oop_returnlists = (
    Vector{Float64}, Vector{Float64},
    ntuple(x -> Vector{dualT}, length(oop_arglists) - 2)...,
)

function wrapfun_oop(ff, inputs::Tuple = ())
    if !isempty(inputs)
        IT = Tuple{map(typeof, inputs)...}
        if IT ∉ NORECOMPILE_OOP_SUPPORTED_ARGS
            throw(NoRecompileArgumentError(IT))
        end
    end
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        ff, oop_arglists,
        oop_returnlists
    )
end

function wrapfun_iip(
        ff,
        inputs::Tuple{T1, T2, T3, T4}
    ) where {T1, T2, T3, T4}
    T = eltype(T2)
    dualT = dualgen(T)
    dualT1 = ArrayInterface.promote_eltype(T1, dualT)
    dualT2 = ArrayInterface.promote_eltype(T2, dualT)
    dualT4 = dualgen(promote_type(T, T4))

    iip_arglists = (
        Tuple{T1, T2, T3, T4},
        Tuple{dualT1, dualT2, T3, T4},
        Tuple{dualT1, T2, T3, dualT4},
        Tuple{dualT1, dualT2, T3, dualT4},
    )

    iip_returnlists = ntuple(x -> Nothing, 4)

    fwt = map(iip_arglists, iip_returnlists) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    return FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

# 3-arg version: compile FunctionWrapper variants with the specified chunk size.
# Uses chunk=CS for u-related duals (Jacobian computation) and chunk=1 for
# t-related duals (time derivative is always scalar, so chunk=1).
function wrapfun_iip(
        ff,
        inputs::Tuple{T1, T2, T3, T4},
        ::Val{CS}
    ) where {T1, T2, T3, T4, CS}
    T = eltype(T2)

    # Jacobian (u-derivative) uses chunk=CS
    dualT_jac = dualgen(T, Val(CS))
    dualT1_jac = ArrayInterface.promote_eltype(T1, dualT_jac)
    dualT2_jac = ArrayInterface.promote_eltype(T2, dualT_jac)

    # Time derivative uses chunk=1 (scalar differentiation w.r.t. t)
    dualT_time = dualgen(T)
    dualT1_time = ArrayInterface.promote_eltype(T1, dualT_time)
    dualT4_time = dualgen(promote_type(T, T4))

    iip_arglists = (
        Tuple{T1, T2, T3, T4},                                  # plain
        Tuple{dualT1_jac, dualT2_jac, T3, T4},                  # Jacobian (u dual, chunk=CS)
        Tuple{dualT1_time, T2, T3, dualT4_time},                # time derivative (chunk=1)
        Tuple{dualT1_jac, dualT2_jac, T3, dualT4_time},         # both
    )

    iip_returnlists = ntuple(x -> Nothing, 4)

    fwt = map(iip_arglists, iip_returnlists) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    return FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

const iip_arglists_default = (
    Tuple{
        Vector{Float64}, Vector{Float64}, Vector{Float64},
        Float64,
    },
    Tuple{
        Vector{Float64}, Vector{Float64},
        SciMLBase.NullParameters,
        Float64,
    },
    Tuple{Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT},
    Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, dualT},
    Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64},
    Tuple{
        Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters,
        Float64,
    },
    Tuple{
        Vector{dualT}, Vector{Float64},
        SciMLBase.NullParameters, dualT,
    },
)
const iip_returnlists_default = ntuple(x -> Nothing, length(iip_arglists_default))

function wrapfun_iip(@nospecialize(ff))
    fwt = map(iip_arglists_default, iip_returnlists_default) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    return FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual}, p, tspan, prob, kwargs)
    if (haskey(kwargs, :callback) && has_continuous_callback(kwargs[:callback])) ||
            (haskey(prob.kwargs, :callback) && has_continuous_callback(prob.kwargs[:callback]))
        return _promote_tspan(eltype(u0).(tspan), kwargs)
    else
        return _promote_tspan(tspan, kwargs)
    end
end

function promote_tspan(
        u0::AbstractArray{<:Complex{<:ForwardDiff.Dual}}, p, tspan, prob,
        kwargs
    )
    return _promote_tspan(real(eltype(u0)).(tspan), kwargs)
end

function promote_tspan(
        u0::AbstractArray{<:ForwardDiff.Dual}, p,
        tspan::Tuple{<:ForwardDiff.Dual, <:ForwardDiff.Dual}, prob, kwargs
    )
    return _promote_tspan(tspan, kwargs)
end
# Copy of the other prob2dtmin dispatch, just for optionality
function prob2dtmin(tspan, ::ForwardDiff.Dual, use_end_time)
    t1, t2 = tspan
    isfinite(t1) || throw(ArgumentError("t0 in the tspan `(t0, t1)` must be finite"))
    if use_end_time && isfinite(t2 - t1)
        return max(eps(t2), eps(t1))
    else
        return max(eps(typeof(t1)), eps(t1))
    end
end

function hasdualpromote(u0, t::Number)
    return hasmethod(
        ArrayInterface.promote_eltype,
        Tuple{Type{typeof(u0)}, Type{dualgen(eltype(u0))}}
    ) &&
        hasmethod(
        promote_rule,
        Tuple{Type{eltype(u0)}, Type{dualgen(eltype(u0))}}
    ) &&
        hasmethod(
        promote_rule,
        Tuple{Type{eltype(u0)}, Type{typeof(t)}}
    )
end

@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::Any) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(
        u::AbstractArray{<:ForwardDiff.Dual{Tag, T}},
        t::Any
    ) where {Tag, T}
    return sqrt(DiffEqBase.__sum(sse, u; init = sse(zero(T))) / DiffEqBase.totallength(u))
end
@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::ForwardDiff.Dual) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(
        u::AbstractArray{<:ForwardDiff.Dual{Tag, T}},
        ::ForwardDiff.Dual
    ) where {Tag, T}
    return sqrt(DiffEqBase.__sum(sse, u; init = sse(zero(T))) / DiffEqBase.totallength(u))
end

if !hasmethod(nextfloat, Tuple{ForwardDiff.Dual})
    # Type piracy. Should upstream
    function Base.nextfloat(d::ForwardDiff.Dual{T, V, N}) where {T, V, N}
        return ForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
    end
    function Base.prevfloat(d::ForwardDiff.Dual{T, V, N}) where {T, V, N}
        return ForwardDiff.Dual{T}(prevfloat(d.value), d.partials)
    end
end

import PrecompileTools
PrecompileTools.@compile_workload begin
    # Scalar operations on Dual numbers (arithmetic, math functions, comparisons)
    d1 = dualT(1.0, ForwardDiff.Partials((0.5,)))
    d2 = dualT(2.0, ForwardDiff.Partials((1.0,)))
    s = 3.14

    # Arithmetic: Dual-Dual and Dual-scalar
    d1 + d2
    d1 - d2
    d1 * d2
    d1 / d2
    d1 + s
    s + d1
    d1 - s
    s - d1
    d1 * s
    s * d1
    d1 / s
    s / d1
    -d1
    abs(d1)

    # Powers and roots
    d1^2
    d1^3
    d2^0.5
    sqrt(d2)
    cbrt(d2)

    # Transcendental functions
    exp(d1)
    log(d2)
    sin(d1)
    cos(d1)
    tan(d1)
    asin(dualT(0.5, ForwardDiff.Partials((1.0,))))
    acos(dualT(0.5, ForwardDiff.Partials((1.0,))))
    atan(d1)
    atan(d1, d2)
    sinh(d1)
    cosh(d1)
    tanh(d1)

    # Comparisons (used in step size control, event detection)
    d1 < d2
    d1 > d2
    d1 <= d2
    d1 >= d2
    d1 == d2
    isnan(d1)
    isinf(d1)
    isfinite(d1)

    # min/max (used in limiters and error control)
    min(d1, d2)
    max(d1, d2)
    min(d1, s)
    max(d1, s)

    # Conversion and promotion
    zero(dualT)
    one(dualT)
    float(d1)
    ForwardDiff.value(d1)
    ForwardDiff.partials(d1)

    # Array operations on Vector{dualT}
    v1 = [d1, d2, dualT(0.0, ForwardDiff.Partials((0.0,)))]
    v2 = [d2, d1, dualT(1.0, ForwardDiff.Partials((0.1,)))]

    # Basic array ops
    v1 + v2
    v1 - v2
    v1 .* v2
    v1 ./ v2
    s .* v1
    v1 .+ s
    v1 .- s
    v1 .^ 2
    v1 .^ 0.5

    # In-place array operations
    out = similar(v1)
    out .= v1 .+ v2
    out .= v1 .- v2
    out .= v1 .* v2
    out .= s .* v1
    out .= v1 .* s .+ v2
    out .= v1 .* s .- v2 .* s

    # Reductions (used in norm calculations, error estimation)
    sum(v1)
    sum(abs2, v1)
    maximum(abs, v1)

    # LinearAlgebra operations
    using LinearAlgebra
    dot(v1, v2)
    norm(v1)
    norm(v1, Inf)
    norm(v1, 1)

    # copy / fill
    copy(v1)
    fill!(out, zero(dualT))

    # SubArray primitive broadcast operations for Float64 and Dual types.
    # These are generic building blocks used by any ODE function with views.
    # Note: fused multi-operand broadcast expressions (e.g. `dy .= k .* y1 .+ k .* y2 .* y3`)
    # create unique nested Broadcasted types per expression and cannot be generically precompiled.
    for T in (Float64, dualT)
        x = zeros(T, 4)
        dx = zeros(T, 4)
        sv1 = @view x[1:2]
        sv2 = @view x[3:4]
        dsv1 = @view dx[1:2]
        k = 0.04

        # Primitive SubArray broadcast operations
        dsv1 .= sv1
        dsv1 .= k .* sv1
        dsv1 .= sv1 .* sv2
        dsv1 .= sv1 .+ sv2
        dsv1 .= sv1 .- sv2
        dsv1 .= sv1 .^ 2
        dsv1 .= .-sv1
    end
end

end
