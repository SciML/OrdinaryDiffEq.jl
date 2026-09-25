module OrdinaryDiffEqExtrapolationPolyesterExt

using Polyester: @batch
import OrdinaryDiffEqExtrapolation: _polyester_foreach

@inline function _polyester_foreach(f, range)
    return @batch for i in range
        f(i)
    end
end

end
