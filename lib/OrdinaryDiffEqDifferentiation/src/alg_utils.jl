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

function DiffEqBase.prepare_alg(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm{0, AD,
                FDT},
            OrdinaryDiffEqImplicitAlgorithm{0, AD, FDT},
            DAEAlgorithm{0, AD, FDT},
            OrdinaryDiffEqExponentialAlgorithm{0, AD, FDT}},
        u0::AbstractArray{T},
        p, prob) where {AD, FDT, T}

    # If not using autodiff or norecompile mode or very large bitsize (like a dual number u0 already)
    # don't use a large chunksize as it will either error or not be beneficial
    if !(alg_autodiff(alg) isa AutoForwardDiff) ||
       (isbitstype(T) && sizeof(T) > 24) ||
       (prob.f isa ODEFunction &&
        prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
        return remake(alg, chunk_size = Val{1}())
    end

    L = StaticArrayInterface.known_length(typeof(u0))
    if L === nothing # dynamic sized
        # If chunksize is zero, pick chunksize right at the start of solve and
        # then do function barrier to infer the full solve
        x = if prob.f.colorvec === nothing
            length(u0)
        else
            maximum(prob.f.colorvec)
        end

        cs = ForwardDiff.pickchunksize(x)
        return remake(alg, chunk_size = Val{cs}())
    else # statically sized
        cs = pick_static_chunksize(Val{L}())
        return remake(alg, chunk_size = cs)
    end
end

@generated function pick_static_chunksize(::Val{chunksize}) where {chunksize}
    x = ForwardDiff.pickchunksize(chunksize)
    :(Val{$x}())
end