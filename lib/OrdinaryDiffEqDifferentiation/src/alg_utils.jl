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

Base.@pure function determine_chunksize(u, alg::SciMLBase.DEAlgorithm)
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
            OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD,
                FDT},
            OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT},
            DAEAlgorithm{CS, AD, FDT},
            OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT}},
        u0::AbstractArray{T},
        p, prob) where {CS, AD, FDT, T}
    prepped_AD = prepare_ADType(alg_autodiff(alg), prob, u0, p, standardtag(alg))

    sparse_prepped_AD = prepare_user_sparsity(prepped_AD, prob)

    # if u0 is a StaticArray or eltype is Complex etc. don't use sparsity
    if (((typeof(u0) <: StaticArray) || (eltype(u0) <: Complex) ||
         (!(prob.f isa DAEFunction) && prob.f.mass_matrix isa MatrixOperator)) &&
        sparse_prepped_AD isa AutoSparse)
        @warn "Input type or problem definition is incompatible with sparse automatic differentiation. Switching to using dense automatic differentiation."
        autodiff = ADTypes.dense_ad(sparse_prepped_AD)
    else
        autodiff = sparse_prepped_AD
    end

    return remake(alg, autodiff = autodiff)
end

function prepare_ADType(autodiff_alg::AutoSparse, prob, u0, p, standardtag)
    SciMLBase.@set autodiff_alg.dense_ad = prepare_ADType(
        ADTypes.dense_ad(autodiff_alg), prob, u0, p, standardtag)
end

function prepare_ADType(autodiff_alg::AutoForwardDiff, prob, u0, p, standardtag)
    tag = if standardtag
        ForwardDiff.Tag(OrdinaryDiffEqTag(), eltype(u0))
    else
        nothing
    end

    T = eltype(u0)

    fwd_cs = OrdinaryDiffEqCore._get_fwd_chunksize_int(autodiff_alg)

    cs = fwd_cs == 0 ? nothing : fwd_cs

    if ((prob.f isa ODEFunction &&
         prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper) ||
        (isbitstype(T) && sizeof(T) > 24)) && (cs == 0 || isnothing(cs))
        return AutoForwardDiff{1}(tag)
    else
        return AutoForwardDiff{cs}(tag)
    end
end

function prepare_ADType(alg::AutoFiniteDiff, prob, u0, p, standardtag)
    # If the autodiff alg is AutoFiniteDiff, prob.f.f isa FunctionWrappersWrapper,
    # and fdtype is complex, fdtype needs to change to something not complex
    if alg.fdtype == Val{:complex}() && (prob.f isa ODEFunction &&
        prob.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
        @warn "AutoFiniteDiff fdtype complex is not compatible with this function"
        return AutoFiniteDiff(fdtype = Val{:forward}())
    end
    return alg
end

function prepare_user_sparsity(ad_alg, prob)
    jac_prototype = prob.f.jac_prototype
    sparsity = prob.f.sparsity

    if !isnothing(sparsity) && !(ad_alg isa AutoSparse)
        if is_sparse_csc(sparsity) && !SciMLBase.has_jac(prob.f)
            if prob.f.mass_matrix isa UniformScaling
                idxs = diagind(sparsity)
                @. @view(sparsity[idxs]) = 1

                if !isnothing(jac_prototype)
                    @. @view(jac_prototype[idxs]) = 1
                end
            else
                idxs = findall(!iszero, prob.f.mass_matrix)
                for idx in idxs
                    sparsity[idx] = prob.f.mass_matrix[idx]
                end

                if !isnothing(jac_prototype)
                    for idx in idxs
                        jac_prototype[idx] = prob.f.mass_matrix[idx]
                    end
                end
            end
        end

        # KnownJacobianSparsityDetector needs an AbstractMatrix
        sparsity = sparsity isa MatrixOperator ? sparsity.A : sparsity

        color_alg = SciMLBase.has_colorvec(prob.f) ?
                    ConstantColoringAlgorithm(
            sparsity, prob.f.colorvec) : GreedyColoringAlgorithm()

        sparsity_detector = ADTypes.KnownJacobianSparsityDetector(sparsity)

        return AutoSparse(
            ad_alg, sparsity_detector = sparsity_detector, coloring_algorithm = color_alg)
    else
        return ad_alg
    end
end

function prepare_ADType(alg::AbstractADType, prob, u0, p, standardtag)
    return alg
end

@generated function pick_static_chunksize(::Val{chunksize}) where {chunksize}
    x = ForwardDiff.pickchunksize(chunksize)
    :(Val{$x}())
end
