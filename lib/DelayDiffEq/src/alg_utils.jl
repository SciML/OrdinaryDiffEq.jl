## SciMLBase Trait Definitions

function SciMLBase.isautodifferentiable(alg::AbstractMethodOfStepsAlgorithm)
    return SciMLBase.isautodifferentiable(alg.alg)
end
function SciMLBase.allows_arbitrary_number_types(alg::AbstractMethodOfStepsAlgorithm)
    return SciMLBase.allows_arbitrary_number_types(alg.alg)
end
function SciMLBase.allowscomplex(alg::AbstractMethodOfStepsAlgorithm)
    return SciMLBase.allowscomplex(alg.alg)
end
SciMLBase.isdiscrete(alg::AbstractMethodOfStepsAlgorithm) = SciMLBase.isdiscrete(alg.alg)
function SciMLBase.isadaptive(alg::AbstractMethodOfStepsAlgorithm)
    inner = alg.alg
    if inner isa Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
        return _sde_isadaptive(inner)
    else
        return SciMLBase.isadaptive(inner)
    end
end

## DelayDiffEq Internal Traits

function isconstrained(alg::AbstractMethodOfStepsAlgorithm{constrained}) where {constrained}
    return constrained
end
OrdinaryDiffEqCore.uses_uprev(alg::AbstractMethodOfStepsAlgorithm, adaptive) = true

function OrdinaryDiffEqCore.isimplicit(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.isimplicit(alg.alg)
end
function OrdinaryDiffEqCore.isdtchangeable(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.isdtchangeable(alg.alg)
end
function OrdinaryDiffEqCore.ismultistep(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.ismultistep(alg.alg)
end
function OrdinaryDiffEqCore.isautoswitch(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.isautoswitch(alg.alg)
end
function OrdinaryDiffEqCore.get_chunksize(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.get_chunksize(alg.alg)
end
function OrdinaryDiffEqCore.get_chunksize_int(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.get_chunksize_int(alg.alg)
end
function OrdinaryDiffEqCore.alg_autodiff(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.alg_autodiff(alg.alg)
end
function OrdinaryDiffEqCore.alg_difftype(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.alg_difftype(alg.alg)
end
function OrdinaryDiffEqCore.standardtag(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.standardtag(alg.alg)
end
function OrdinaryDiffEqCore.concrete_jac(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.concrete_jac(alg.alg)
end
function OrdinaryDiffEqCore.alg_extrapolates(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.alg_extrapolates(alg.alg)
end
function SciMLBase.alg_order(alg::AbstractMethodOfStepsAlgorithm)
    inner = alg.alg
    if inner isa Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
        return _sde_alg_order(inner)
    else
        return SciMLBase.alg_order(inner)
    end
end
function OrdinaryDiffEqCore.alg_adaptive_order(alg::AbstractMethodOfStepsAlgorithm)
    return OrdinaryDiffEqCore.alg_adaptive_order(alg.alg)
end

"""
    iscomposite(alg)

Return if algorithm `alg` is a composite algorithm.
"""
iscomposite(alg) = false
iscomposite(::OrdinaryDiffEqCore.OrdinaryDiffEqCompositeAlgorithm) = true
iscomposite(alg::AbstractMethodOfStepsAlgorithm) = iscomposite(alg.alg)

# SDE algorithms define alg_order in StochasticDiffEqCore's namespace, not SciMLBase's.
# Bridge the two so that alg_maximum_order works for SDE algorithms.
function OrdinaryDiffEqCore.alg_maximum_order(alg::SDEAlgUnion)
    return max(1, floor(Int, _sde_alg_order(alg)))
end
function OrdinaryDiffEqCore.alg_maximum_order(alg::AbstractMethodOfStepsAlgorithm)
    inner = alg.alg
    if inner isa Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
        # SDE alg_order returns Rational; convert safely via floor
        return max(1, floor(Int, _sde_alg_order(inner)))
    else
        return OrdinaryDiffEqCore.alg_maximum_order(inner)
    end
end

function DiffEqBase.prepare_alg(alg::MethodOfSteps, u0, p, prob)
    return MethodOfSteps(
        DiffEqBase.prepare_alg(alg.alg, u0, p, prob);
        constrained = isconstrained(alg),
        fpsolve = alg.fpsolve
    )
end
