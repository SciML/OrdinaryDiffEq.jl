################################################################################
# Ihat2 methods for weak-order cache types
# These extend StochasticDiffEqCore.Ihat2 with cache-type-specific dispatch.

function Ihat2(cache::Union{DRI1ConstantCache, DRI1Cache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return (_dW[k] * _dW[l] - sqdt * _dZ[k]) / 2
    elseif l < k
        return (_dW[k] * _dW[l] + sqdt * _dZ[l]) / 2
    else
        # l == k
        return (_dW[k]^2 - sqdt^2) / 2
    end
end

function Ihat2(cache::Union{RSConstantCache, RSCache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return -_dW[k] * _dZ[l]
    elseif l < k
        return _dW[k] * _dZ[l]
    else
        # l == k
        return zero(_dW[k])
    end
end

function Ihat2(cache::Union{PL1WMConstantCache, PL1WMCache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return -_dZ[Int(1 + 1 // 2 * (l - 3) * l + k)]
    elseif l < k
        return _dZ[Int(1 + 1 // 2 * (k - 3) * k + l)]
    else
        # l == k
        return -sqdt^2
    end
end

function Ihat2(cache::Union{NONConstantCache, NONCache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return _dZ[k]
    elseif l < k
        return _dW[l]
    else
        # l == k
        return _dW[k]
    end
end
