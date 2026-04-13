@muladd function calc_twopoint_random!(_dW, sqdt, dW)
    @.. _dW = ifelse(sign(dW) > false, sqdt, -sqdt)
    return nothing
end

@muladd function calc_twopoint_random(sqdt, dW)
    return ifelse(sign(dW) > false, sqdt, -sqdt)
end

function calc_threepoint_random!(_dW, sq3dt, quantile, dW_scaled)
    @. _dW = ifelse(abs(dW_scaled) > -quantile, ifelse(dW_scaled < quantile, -sq3dt, sq3dt), zero(sq3dt))
    return nothing
end

@muladd function calc_threepoint_random(sq3dt, quantile, dW_scaled)
    return ifelse(abs(dW_scaled) > -quantile, ifelse(dW_scaled < quantile, -sq3dt, sq3dt), zero(sq3dt))
end

# Ihat2 function stub: methods are defined in StochasticDiffEqWeak alongside their cache types.
function Ihat2 end
