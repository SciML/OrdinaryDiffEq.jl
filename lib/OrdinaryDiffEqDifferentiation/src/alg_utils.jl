# Extract AD type parameter from algorithm, returning as Val to ensure type stability for boolean options.
function _alg_autodiff(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have an autodifferentiation option defined.")
end
_alg_autodiff(::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(::DAEAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(::OrdinaryDiffEqImplicitAlgorithm{CS, AD}) where {CS, AD} = Val{AD}()
_alg_autodiff(alg::CompositeAlgorithm) = _alg_autodiff(alg.algs[end])
function _alg_autodiff(::Union{OrdinaryDiffEqExponentialAlgorithm{CS, AD},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD}
}) where {
        CS, AD
}
    Val{AD}()
end

function alg_autodiff(alg)
    autodiff = _alg_autodiff(alg)
    if autodiff == Val(false)
        return AutoFiniteDiff()
    elseif autodiff == Val(true)
        return AutoForwardDiff()
    else
        return _unwrap_val(autodiff)
    end
end

Base.@pure function determine_chunksize(u, alg::DiffEqBase.DEAlgorithm)
    determine_chunksize(u, get_chunksize(alg))
end
Base.@pure function determine_chunksize(u, CS)
    if CS != 0
        return CS
    else
        return ForwardDiff.pickchunksize(length(u))
    end
end
