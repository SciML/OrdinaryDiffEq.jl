"""
    calc_twopoint_random!(_dW, sqdt, dW) -> nothing

In-place form of [`calc_twopoint_random`](@ref), writing the discrete increments into
`_dW`.
"""
@muladd function calc_twopoint_random!(_dW, sqdt, dW)
    @.. _dW = ifelse(sign(dW) > false, sqdt, -sqdt)
    return nothing
end

"""
    calc_twopoint_random(sqdt, dW)

Two-point discrete random variable replacing the Brownian increment in a weak
(moment-matching) scheme.

Returns `±sqdt`, taking the sign of `dW`, so the result has the same mean and
variance as the Brownian increment while taking only two values. Weak-order methods
may substitute such a variable because they only need the moments of the noise to
match, and sampling two points is cheaper than sampling a Gaussian.

`sqdt` is `sqrt(dt)`. See [`calc_threepoint_random`](@ref) for the variant that also
matches the fourth moment.
"""
@muladd function calc_twopoint_random(sqdt, dW)
    return ifelse(sign(dW) > false, sqdt, -sqdt)
end

"""
    calc_threepoint_random!(_dW, sq3dt, quantile, dW_scaled) -> nothing

In-place form of [`calc_threepoint_random`](@ref), writing the discrete increments
into `_dW`.
"""
function calc_threepoint_random!(_dW, sq3dt, quantile, dW_scaled)
    @. _dW = ifelse(abs(dW_scaled) > -quantile, ifelse(dW_scaled < quantile, -sq3dt, sq3dt), zero(sq3dt))
    return nothing
end

"""
    calc_threepoint_random(sq3dt, quantile, dW_scaled)

Three-point discrete random variable replacing the Brownian increment in a weak
(moment-matching) scheme.

Returns `-sq3dt`, `0`, or `+sq3dt` depending on where the standardized increment
`dW_scaled` falls relative to `quantile`. Unlike the two-point variable of
[`calc_twopoint_random`](@ref) this matches the fourth moment of the Gaussian as
well, which weak order 2 schemes require.

`sq3dt` is `sqrt(3dt)`, and `quantile` is the (negative) standard-normal quantile
that gives the zero outcome its correct probability.
"""
@muladd function calc_threepoint_random(sq3dt, quantile, dW_scaled)
    return ifelse(abs(dW_scaled) > -quantile, ifelse(dW_scaled < quantile, -sq3dt, sq3dt), zero(sq3dt))
end

"""
    Ihat2(...)

Approximation of the second-order multiple stochastic integral used by the weak
order 2 schemes.

Only the generic function lives here; the methods are defined in
StochasticDiffEqWeak alongside the cache types they dispatch on, because their
argument lists differ between the DRI/RI (Rößler) and RDI families.
"""
function Ihat2 end
