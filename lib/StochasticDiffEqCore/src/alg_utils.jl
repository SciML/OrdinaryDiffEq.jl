## SciMLBase Trait Definitions

SciMLBase.allows_late_binding_tstops(::StochasticDiffEqAlgorithm) = true
SciMLBase.allows_late_binding_tstops(::StochasticDiffEqRODEAlgorithm) = true

SciMLBase.supports_solve_rng(
    ::SciMLBase.AbstractSDEProblem,
    ::StochasticDiffEqAlgorithm,
) = true

SciMLBase.supports_solve_rng(
    ::SciMLBase.AbstractRODEProblem,
    ::StochasticDiffEqRODEAlgorithm,
) = true

function SciMLBase.isautodifferentiable(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return true
end
function SciMLBase.allows_arbitrary_number_types(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return true
end
function SciMLBase.allowscomplex(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return true
end
SciMLBase.isdiscrete(alg::StochasticDiffEqJumpAlgorithm) = true

function SciMLBase.forwarddiffs_model(
        alg::Union{
            StochasticDiffEqNewtonAlgorithm,
            StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
            StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
        }
    )
    return OrdinaryDiffEqCore.alg_autodiff(alg) isa ADTypes.AutoForwardDiff
end

# Required for initialization, because ODECore._initialize_dae! calls it during
# OverrideInit
function OrdinaryDiffEqCore.has_autodiff(
        ::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
            StochasticDiffEqJumpAlgorithm,
        }
    )
    return false
end
for T in [
        StochasticDiffEqNewtonAlgorithm, StochasticDiffEqNewtonAdaptiveAlgorithm,
        StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
        StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
    ]
    @eval OrdinaryDiffEqCore.has_autodiff(::$T) = true
end

_alg_autodiff(::StochasticDiffEqNewtonAlgorithm{T, AD}) where {T, AD} = Val{AD}()
_alg_autodiff(::StochasticDiffEqNewtonAdaptiveAlgorithm{T, AD}) where {T, AD} = Val{AD}()
function _alg_autodiff(::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{T, AD}) where {T, AD}
    return Val{AD}()
end
function _alg_autodiff(
        ::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            T, AD,
        }
    ) where {T, AD}
    return Val{AD}()
end
_alg_autodiff(alg::StochasticCompositeAlgorithm) = _alg_autodiff(alg.algs[end])

function OrdinaryDiffEqCore.alg_autodiff(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    ad = _alg_autodiff(alg)
    if ad == Val(false)
        return ADTypes.AutoFiniteDiff()
    elseif ad == Val(true)
        return ADTypes.AutoForwardDiff()
    else
        return SciMLBase._unwrap_val(ad)
    end
end

isadaptive(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false
function isadaptive(
        alg::Union{
            StochasticDiffEqAdaptiveAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
            StochasticDiffEqJumpAdaptiveAlgorithm,
            StochasticDiffEqJumpDiffusionAdaptiveAlgorithm,
        }
    )
    return true
end
function isadaptive(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return all(isadaptive.(alg.algs))
end
function isadaptive(
        prob, alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return isadaptive(alg)
end

## StochasticDiffEq Internal Traits

function qmax_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm})
    return isadaptive(alg) ? 9 // 8 : 0
end
function qmin_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm})
    return isadaptive(alg) ? 1 // 5 : 0
end

delta_default(alg) = 1 // 1

function OrdinaryDiffEqCore.gamma_default(
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
    )
    return isadaptive(alg) ? 9 // 10 : 0
end

ispredictive(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false
isstandard(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false
function qsteady_min_default(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return 1
end
function qsteady_max_default(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return 1
end

# Extend ODE's isaposteriori for SDE algorithms — default is false.
# Solver subpackages may extend for specific algorithms (e.g., CaoTauLeaping).

# SDE algorithms never extrapolate (no uprev2 storage).
OrdinaryDiffEqCore.alg_extrapolates(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false

# SDE algorithms are never FSAL (update_fsal! is a no-op for SDE).
OrdinaryDiffEqCore.isfsal(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = false

# Composite alg_order fallback
function alg_order(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return maximum(alg_order.(alg.algs))
end
get_current_alg_order(alg::StochasticDiffEqAlgorithm, cache) = alg_order(alg)
get_current_alg_order(alg::StochasticDiffEqRODEAlgorithm, cache) = alg_order(alg)
function get_current_alg_order(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        },
        cache
    )
    return alg_order(alg.algs[cache.current])
end

function beta2_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm})
    return isadaptive(alg) ? 2 // (5alg_order(alg)) : 0
end
function beta1_default(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}, beta2)
    return isadaptive(alg) ? 7 // (10alg_order(alg)) : 0
end

isdtchangeable(alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}) = true

# Default alg_interpretation — Ito for SDE.
# Solver subpackages override for Stratonovich methods.
function SciMLBase.alg_interpretation(alg::StochasticDiffEqAlgorithm)
    return SciMLBase.AlgorithmInterpretation.Ito
end

# Default alg_compatible — false for abstract SDE type, true for RODE
function alg_compatible(
        prob, alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return true
end
alg_compatible(prob, alg::StochasticDiffEqAlgorithm) = false

function alg_compatible(
        prob::DiffEqBase.AbstractSDEProblem,
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return max((alg_compatible(prob, a) for a in alg.algs)...)
end

# Trait: whether an algorithm supports regular jumps in a JumpProblem.
# Default is false; EM and ImplicitEM override to true in their subpackages.
supports_regular_jumps(alg) = false

# JumpProblem compatibility defaults
function alg_compatible(prob::JumpProblem, alg::StochasticDiffEqAlgorithm)
    return alg_compatible(prob.prob, alg) &&
        (supports_regular_jumps(alg) || prob.regular_jump === nothing) &&
        prob.prob isa DiffEqBase.AbstractSDEProblem
end

function alg_compatible(
        prob::JumpProblem,
        alg::Union{StochasticDiffEqJumpAdaptiveAlgorithm, StochasticDiffEqJumpAlgorithm}
    )
    return prob.prob isa DiscreteProblem
end

function alg_needs_extra_process(
        alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return false
end
function alg_needs_extra_process(
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return max((alg_needs_extra_process(a) for a in alg.algs)...)
end

function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqNewtonAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end
function OrdinaryDiffEqDifferentiation._alg_autodiff(
        alg::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val{AD}()
end

function OrdinaryDiffEqCore.get_current_alg_autodiff(alg::StochasticDiffEqCompositeAlgorithm, cache)
    return OrdinaryDiffEqCore.alg_autodiff(alg.algs[cache.current])
end

function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqNewtonAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end
function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end
function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end
function OrdinaryDiffEqCore.get_chunksize(
        alg::StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
            CS, AD, FDT, ST, CJ, Controller,
        }
    ) where {CS, AD, FDT, ST, CJ, Controller}
    return Val(CS)
end

@static if isdefined(OrdinaryDiffEqCore, :standardtag)
    OrdinaryDiffEqCore.standardtag(
        alg::Union{
            StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller},
            StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller},
        }
    ) where {CS, AD, FDT, ST, CJ, Controller} = ST
end

@static if isdefined(OrdinaryDiffEqCore, :alg_difftype)
    OrdinaryDiffEqCore.alg_difftype(
        alg::Union{
            StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller},
            StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller},
        }
    ) where {
        CS, AD, FDT, ST, CJ, Controller,
    } = FDT
end

@static if isdefined(OrdinaryDiffEqCore, :concrete_jac)
    OrdinaryDiffEqCore.concrete_jac(
        alg::Union{
            StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller},
            StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller},
        }
    ) where {
        CS, AD, FDT, ST, CJ, Controller,
    } = CJ
end

alg_mass_matrix_compatible(alg::StochasticDiffEqAlgorithm) = false
alg_mass_matrix_compatible(alg::StochasticDiffEqRODEAlgorithm) = false
alg_can_repeat_jac(alg::StochasticDiffEqAlgorithm) = true

function alg_mass_matrix_compatible(
        alg::Union{
            StochasticDiffEqNewtonAlgorithm, StochasticDiffEqNewtonAdaptiveAlgorithm,
        }
    )
    if alg.symplectic
        return true
    elseif alg.theta == 1
        return true
    else
        error("Algorithm must be set as symplectic or theta=1 for mass matrices")
    end
end

is_split_step(::StochasticDiffEqAlgorithm) = false

alg_stability_size(alg::StochasticDiffEqAlgorithm) = 0 # default, overridden per-alg

# is_composite_algorithm trait is defined in OrdinaryDiffEqCore and extended in
# integrator_utils.jl for SDE composite algorithm types.
function unwrap_alg(integrator, is_nlsolve)
    alg = integrator.alg
    if !is_composite_algorithm(alg)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        num = is_nlsolve ? 2 : 1
        if num == 1
            return alg.algs[1]
        elseif num == 2
            return alg.algs[2]
        else
            return alg.algs[num]
        end
    else
        if integrator.cache.current == 1
            return alg.algs[1]
        elseif integrator.cache.current == 2
            return alg.algs[2]
        else
            return alg.algs[integrator.cache.current]
        end
    end
end

alg_control_rate(::StochasticDiffEqAlgorithm) = false
alg_control_rate(::StochasticDiffEqRODEAlgorithm) = false
