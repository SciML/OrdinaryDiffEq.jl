# Extract AD type parameter from algorithm, returning as Val to ensure type stability for boolean options.
function _alg_autodiff(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have an autodifferentiation option defined.")
end
function _alg_autodiff(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD}
    alg.autodiff
end
_alg_autodiff(alg::DAEAlgorithm{CS, AD}) where {CS, AD} = alg.autodiff
_alg_autodiff(alg::OrdinaryDiffEqImplicitAlgorithm{CS, AD}) where {CS, AD} = alg.autodiff
_alg_autodiff(alg::CompositeAlgorithm) = _alg_autodiff(alg.algs[end])
function _alg_autodiff(alg::Union{OrdinaryDiffEqExponentialAlgorithm{CS, AD},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD}
}) where {
        CS, AD
}
    alg.autodiff
end

function alg_autodiff(alg)
    autodiff = _alg_autodiff(alg)

    if autodiff == Val(true)
        return AutoForwardDiff()
    elseif autodiff == Val(false)
        return AutoFiniteDiff()
    else
        return autodiff
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
    # If prob.f.f is a FunctionWrappersWrappers from ODEFunction, need to set chunksize to 1

    if alg_autodiff(alg) isa AutoForwardDiff && ((prob.f isa ODEFunction &&
         prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper) ||
        (isbitstype(T) && sizeof(T) > 24))
        return remake(
            alg, autodiff = AutoForwardDiff(chunksize = 1, tag = alg_autodiff(alg).tag))
    end

    # If the autodiff alg is AutoFiniteDiff, prob.f.f isa FunctionWrappersWrapper,
    # and fdtype is complex, fdtype needs to change to something not complex
    if alg_autodiff(alg) isa AutoFiniteDiff
        if alg_difftype(alg) == Val{:complex} && (prob.f isa ODEFunction &&
            prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
            @warn "AutoFiniteDiff fdtype complex is not compatible with this function"
            return remake(alg, autodiff = AutoFiniteDiff(fdtype = Val{:forward}()))
        end
        return alg
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
        return remake(alg,
            autodiff = AutoForwardDiff(
                chunksize = cs))
    else # statically sized
        cs = pick_static_chunksize(Val{L}())
        cs = SciMLBase._unwrap_val(cs)
        return remake(
            alg, autodiff = AutoForwardDiff(chunksize = cs))
    end
end

@generated function pick_static_chunksize(::Val{chunksize}) where {chunksize}
    x = ForwardDiff.pickchunksize(chunksize)
    :(Val{$x}())
end
