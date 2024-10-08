const DEFAULT_AUTODIFF = AutoForwardDiff()

# Extract AD type parameter from algorithm, returning as Val to ensure type stability for boolean options.
function alg_autodiff(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have an autodifferentiation option defined.")
end
alg_autodiff(::OrdinaryDiffEqAdaptiveImplicitAlgorithm{AD}) where {AD} = AD
alg_autodiff(::DAEAlgorithm{AD}) where {AD} = AD
alg_autodiff(::OrdinaryDiffEqImplicitAlgorithm{AD}) where {AD} = AD
alg_autodiff(alg::CompositeAlgorithm) = _alg_autodiff(alg.algs[end])
function alg_autodiff(::Union{OrdinaryDiffEqExponentialAlgorithm{AD},
        OrdinaryDiffEqAdaptiveExponentialAlgorithm{AD}
}) where {AD}
    AD
end

set_chunksize(backend::AbstractADType, ::Val{C}) where {C} = backend

function set_chunksize(backend::AutoForwardDiff{<:Any, T}, ::Val{C}) where {C, T}
    return AutoForwardDiff{C, T}(backend.tag)
end

clever_chunksize(backend::AbstractADType, x::AbstractArray) = backend

function DiffEqBase.prepare_alg(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm{AD, FDT},
            OrdinaryDiffEqImplicitAlgorithm{AD, FDT},
            DAEAlgorithm{AD, FDT},
            OrdinaryDiffEqExponentialAlgorithm{AD, FDT}},
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
