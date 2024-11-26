# Extract AD type parameter from algorithm, returning as Val to ensure type stability for boolean options.
function _alg_autodiff(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have an autodifferentiation option defined.")
end
_alg_autodiff(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD}) where {CS, AD} = alg.autodiff
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


    autodiff = prepare_ADType(alg_autodiff(alg), prob, u0, p, standardtag(alg))

    #sparsity preparation

    sparsity = prob.f.sparsity

    if sparsity isa SparseMatrixCSC
        if f.mass_matrix isa UniformScaling
            idxs = diagind(sparsity)
            @. @view(sparsity[idxs]) = 1
        else
            idxs = findall(!iszero, f.mass_matrix)
            @. @view(sparsity[idxs]) = @view(f.mass_matrix[idxs])
        end
    end

    sparsity_detector = isnothing(sparsity) ? TracerSparsityDetector() : ADTypes.KnownJacobianSparsityDetector(sparsity)
    color_alg = DiffEqBase.has_colorvec(prob.f) ? ADTypes.ConstantColoringAlgorithm(sparsity, prob.f.colorvec) : GreedyColoringAlgorithm()

    autodiff = AutoSparse(autodiff, sparsity_detector = sparsity_detector, coloring_algorithm = color_alg)

    alg = remake(alg, autodiff = autodiff)

    return alg
end

function prepare_ADType(autodiff_alg::AutoSparse, prob, u0, p, standardtag)
    prepare_ADType(dense_ad(autodiff_alg), prob, u0, p, standardtag)
end

function prepare_ADType(autodiff_alg::AutoForwardDiff, prob, u0, p, standardtag)
    tag = if standardtag
        ForwardDiff.Tag(OrdinaryDiffEqTag(), eltype(u0))
    else
        nothing
    end

    T = eltype(u0)

    if ((prob.f isa ODEFunction &&
      prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper) ||
     (isbitstype(T) && sizeof(T) > 24))
        autodiff_alg = AutoForwardDiff(chunksize = 1, tag = tag)
    end

    #L = StaticArrayInterface.known_length(typeof(u0))
    #if L === nothing # dynamic sized
        # If chunksize is zero, pick chunksize right at the start of solve and
        # then do function barrier to infer the full solve
    #    x = if prob.f.colorvec === nothing
    #        length(u0)
    #    else
    #        maximum(prob.f.colorvec)
    #    end

    #    cs = ForwardDiff.pickchunksize(x)
    #    return remake(alg,
    #        autodiff = AutoForwardDiff(
    #            chunksize = cs, tag = tag))
    #else # statically sized
    #    cs = pick_static_chunksize(Val{L}())
    #    cs = SciMLBase._unwrap_val(cs)
    #    return remake(
    #        alg, autodiff = AutoForwardDiff(chunksize = cs, tag = tag))
    #end
    autodiff_alg
end

function prepare_ADType(alg::AutoFiniteDiff, prob, u0, p, standardtag)
    # If the autodiff alg is AutoFiniteDiff, prob.f.f isa FunctionWrappersWrapper,
    # and fdtype is complex, fdtype needs to change to something not complex
    if alg.fdtype == Val{:complex}() && (prob.f isa ODEFunction && prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
         @warn "AutoFiniteDiff fdtype complex is not compatible with this function"
         return AutoFiniteDiff(fdtype = Val{:forward}())
    end
    return alg
end

function prepare_ADType(alg::AbstractADType, prob, u0,p,standardtag)
    return alg
end

#function prepare_ADType(alg::DiffEqAutoAD, prob, u0, p, standardtag)

#end

@generated function pick_static_chunksize(::Val{chunksize}) where {chunksize}
    x = ForwardDiff.pickchunksize(chunksize)
    :(Val{$x}())
end
