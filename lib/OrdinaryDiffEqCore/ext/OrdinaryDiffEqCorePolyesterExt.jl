module OrdinaryDiffEqCorePolyesterExt

using Polyester: @batch
import OrdinaryDiffEqCore: _polyester_batch

@inline function _polyester_batch(f, range)
    @batch for i in range
        f(i)
    end
end

end
