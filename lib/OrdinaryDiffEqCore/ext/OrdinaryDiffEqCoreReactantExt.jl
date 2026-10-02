module OrdinaryDiffEqCoreReactantExt

using OrdinaryDiffEqCore: OrdinaryDiffEqCore
using Reactant: Reactant, TracedRNumber

bits(x::TracedRNumber{T}) where {T} = Reactant.Ops.bitcast_convert(Base.uinttype(T), x)

magnitude_bits(x::TracedRNumber{T}) where {T} = bits(x) & ~Base.sign_mask(T)

# Compiled code flushes subnormal floats to zero, so the comparison is made on bit patterns.
# The bits of `eps(t)` for finite `t` are `2^(E - 1)` for a biased exponent `E <= p`, where
# `eps(t)` is subnormal, and `(E - p) << p` otherwise. Non-negative floats and non-NaN bit
# patterns order alike, and a NaN `dt` compares above every finite `eps`.
function OrdinaryDiffEqCore.dt_below_time_eps(
        dt::Union{T, TracedRNumber{T}}, t::TracedRNumber{T}
    ) where {T <: Base.IEEEFloat}
    U = Base.uinttype(T)
    p = U(Base.significand_bits(T))
    biased = div(bits(t) & Base.exponent_mask(T), one(U) << p)
    subnormal_eps = Reactant.Ops.shift_left(
        Reactant.promote_to(TracedRNumber{U}, one(U)), max(biased, one(U)) - one(U)
    )
    eps_bits = ifelse(biased > p, (biased - p) * (one(U) << p), subnormal_eps)
    dt_bits = magnitude_bits(Reactant.promote_to(TracedRNumber{T}, dt))
    return (biased != Base.exponent_mask(T) >> p) & (dt_bits <= eps_bits)
end

# `nextfloat(a) - a` is exact, `eps(floatmax(T))` caps `a == floatmax(T)`, and a
# non-finite `x` gives NaN. A subnormal spacing, which compiled code would flush to zero,
# is raised to `floatmin(T)`, the nearest step that compiled code can take.
function OrdinaryDiffEqCore.value_eps(x::TracedRNumber{T}) where {T <: Base.IEEEFloat}
    a = abs(x)
    return max(min(nextfloat(a) - a, eps(floatmax(T))), floatmin(T))
end

# `FastPower.fastlog2`: the significand and exponent are read from the bits of `x`.
function fastlog2(x::TracedRNumber{Float32})
    a = 0.338953f0
    b = 1.8596461f0
    c = 0.523692f0
    ux1i = bits(x)
    upper = iszero(ux1i & 0x00400000)
    biased = Reactant.Ops.convert(
        TracedRNumber{Float32}, div(ux1i & 0x7F800000, 0x00800000)
    )
    exponent = biased - ifelse(upper, 127.0f0, 126.0f0)
    ux2i = (ux1i & 0x007FFFFF) | ifelse(upper, 0x3f800000, 0x3f000000)
    signif = Reactant.Ops.bitcast_convert(Float32, ux2i)
    quot = (signif * a + b) / (signif + c)
    return (signif - 1.0f0) * quot + exponent
end

# `FastPower.fastpower(x::T, y::T)`, which the controllers use for floats on the host, up to
# the rounding of `exp2`. An `x` whose `Float32` value is subnormal is flushed to zero.
function fastpower(x::TracedRNumber{T}, y::TracedRNumber{T}) where {T}
    x32 = Reactant.Ops.convert(TracedRNumber{Float32}, x)
    y32 = Reactant.Ops.convert(TracedRNumber{Float32}, y)
    approx = Reactant.Ops.convert(TracedRNumber{T}, 2.0f0^(y32 * fastlog2(x32)))
    return ifelse(
        iszero(x), zero(approx),
        ifelse(isinf(x) & isinf(y), T(Inf), approx)
    )
end

promote_traced(x, ::Type{T}) where {T} = Reactant.promote_to(TracedRNumber{T}, x)

OrdinaryDiffEqCore.controller_fastpower(
    x::TracedRNumber{T}, y::Union{T, TracedRNumber{T}}
) where {T <: Base.IEEEFloat} = fastpower(x, promote_traced(y, T))
OrdinaryDiffEqCore.controller_fastpower(
    x::T, y::TracedRNumber{T}
) where {T <: Base.IEEEFloat} = fastpower(promote_traced(x, T), y)

end
