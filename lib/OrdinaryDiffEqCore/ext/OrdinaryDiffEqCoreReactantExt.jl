module OrdinaryDiffEqCoreReactantExt

using OrdinaryDiffEqCore: OrdinaryDiffEqCore
using ReactantCore: promote_to_traced
# Neither name has a public equivalent: methods for traced numbers dispatch on
# `TracedRNumber`, and Reactant defines no `reinterpret` for traced numbers.
using Reactant: TracedRNumber
using Reactant.Ops: bitcast_convert

# Compiled code flushes subnormal floats to zero in arithmetic and comparisons, but not in
# bit patterns, so the functions below work on bits wherever the host result depends on a
# subnormal value.

const IEEEFloat = Union{Float16, Float32, Float64}

uint(::Type{Float16}) = UInt16
uint(::Type{Float32}) = UInt32
uint(::Type{Float64}) = UInt64

significand_bits(::Type{T}) where {T} = uint(T)(precision(T) - 1)

magnitude_mask(::Type{T}) where {T} = typemax(uint(T)) >> 1

exponent_mask(::Type{T}) where {T} =
    magnitude_mask(T) & ~((one(uint(T)) << significand_bits(T)) - one(uint(T)))

bits(x::TracedRNumber{T}) where {T} = bitcast_convert(uint(T), x)

magnitude_bits(x::TracedRNumber{T}) where {T} = bits(x) & magnitude_mask(T)

# The bits of `eps(t)` for finite `t` are `2^(E - 1)` for a biased exponent `E <= p`, where
# `eps(t)` is subnormal, and `(E - p) << p` otherwise. Non-negative floats and non-NaN bit
# patterns order alike, and a NaN `dt` compares above every finite `eps`.
function OrdinaryDiffEqCore.dt_below_time_eps(
        dt::Union{T, TracedRNumber{T}}, t::TracedRNumber{T}
    ) where {T <: IEEEFloat}
    U = uint(T)
    p = significand_bits(T)
    biased = div(bits(t) & exponent_mask(T), one(U) << p)
    subnormal_eps = promote_to_traced(U(2))^(max(biased, one(U)) - one(U))
    eps_bits = ifelse(biased > p, (biased - p) * (one(U) << p), subnormal_eps)
    dt_bits = magnitude_bits(promote_to_traced(dt))
    return (biased != exponent_mask(T) >> p) & (dt_bits <= eps_bits)
end

# `nextfloat(a) - a` is exact, `eps(floatmax(T))` caps `a == floatmax(T)`, and a
# non-finite `x` gives NaN. A subnormal spacing, which compiled code would flush to zero,
# is raised to `floatmin(T)`, the nearest step that compiled code can take.
function OrdinaryDiffEqCore.value_eps(x::TracedRNumber{T}) where {T <: IEEEFloat}
    a = abs(x)
    return max(min(nextfloat(a) - a, eps(floatmax(T))), floatmin(T))
end

# The bits of `abs(convert(Float32, x))` as the host computes them. Below `floatmin(Float32)`
# the conversion is subnormal, which compiled conversion flushes, so the significand is
# rounded to the subnormal grid (nearest, ties to even) directly.
float32_magnitude_bits(x::TracedRNumber{Float32}) = magnitude_bits(x)
function float32_magnitude_bits(x::TracedRNumber{Float64})
    a = abs(x)
    normal = magnitude_bits(convert(TracedRNumber{Float32}, a))
    subnormal = unsafe_trunc(UInt32, round(a * 0x1p149))
    return ifelse(a < floatmin(Float32), subnormal, normal)
end

# `FastPower.fastlog2` of the `Float32` with magnitude bits `ux1i`.
function fastlog2(ux1i::TracedRNumber{UInt32})
    a = 0.338953f0
    b = 1.8596461f0
    c = 0.523692f0
    upper = iszero(ux1i & 0x00400000)
    biased = convert(TracedRNumber{Float32}, div(ux1i & 0x7F800000, 0x00800000))
    exponent = biased - ifelse(upper, 127.0f0, 126.0f0)
    ux2i = (ux1i & 0x007FFFFF) | ifelse(upper, 0x3f800000, 0x3f000000)
    signif = bitcast_convert(Float32, ux2i)
    quot = (signif * a + b) / (signif + c)
    return (signif - 1.0f0) * quot + exponent
end

# `FastPower.fastpower(x::T, y::T)`, which the controllers use for floats on the host, up to
# the rounding of `exp2` (within `eps(Float32)` relative). The host's `exp2` returns zero
# below `2^-126`, as the compiled one does after flushing, so results are never subnormal.
# Zero and infinity are classified on bits too: the compiler otherwise rewrites the tests
# into float comparisons, which treat a subnormal `x` as zero.
function fastpower(x::TracedRNumber{T}, y::TracedRNumber{T}) where {T}
    y32 = convert(TracedRNumber{Float32}, y)
    approx = convert(TracedRNumber{T}, 2.0f0^(y32 * fastlog2(float32_magnitude_bits(x))))
    xbits = magnitude_bits(x)
    both_inf = (xbits == exponent_mask(T)) & (magnitude_bits(y) == exponent_mask(T))
    return ifelse(iszero(xbits), zero(approx), ifelse(both_inf, T(Inf), approx))
end

const ControllerFloat = Union{Float32, Float64}

OrdinaryDiffEqCore.controller_fastpower(
    x::TracedRNumber{T}, y::Union{T, TracedRNumber{T}}
) where {T <: ControllerFloat} = fastpower(x, promote_to_traced(y))
OrdinaryDiffEqCore.controller_fastpower(
    x::T, y::TracedRNumber{T}
) where {T <: ControllerFloat} = fastpower(promote_to_traced(x), y)

unsupported_float16() = throw(
    ArgumentError("Float16 step-size control is not supported inside Reactant compilation")
)
OrdinaryDiffEqCore.controller_fastpower(
    x::TracedRNumber{Float16}, y::Union{Float16, TracedRNumber{Float16}}
) = unsupported_float16()
OrdinaryDiffEqCore.controller_fastpower(x::Float16, y::TracedRNumber{Float16}) =
    unsupported_float16()

end
