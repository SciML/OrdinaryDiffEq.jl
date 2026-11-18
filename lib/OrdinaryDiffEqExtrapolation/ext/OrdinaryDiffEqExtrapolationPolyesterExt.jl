module OrdinaryDiffEqExtrapolationPolyesterExt

using Polyester: @batch
import OrdinaryDiffEqExtrapolation: _polyester_foreach, _nchunks, _chunk

@inline function _polyester_foreach(f, items)
    nchunks = _nchunks(items)
    return @batch for c in 1:nchunks
        for i in _chunk(items, c, nchunks)
            f(i, c)
        end
    end
end

end
