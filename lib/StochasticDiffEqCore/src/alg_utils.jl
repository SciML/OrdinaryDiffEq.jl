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


function OrdinaryDiffEqCore.alg_autodiff(
        alg::Union{
            StochasticDiffEqNewtonAlgorithm,
            StochasticDiffEqNewtonAdaptiveAlgorithm,
            StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
            StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
        }
    )
    ad = alg.autodiff
    if ad == Val(false)
        return ADTypes.AutoFiniteDiff()
    elseif ad == Val(true)
        return ADTypes.AutoForwardDiff()
    else
        return SciMLBase._unwrap_val(ad)
    end
end
OrdinaryDiffEqCore.alg_autodiff(alg::StochasticDiffEqCompositeAlgorithm) = OrdinaryDiffEqCore.alg_autodiff(alg.algs[end])

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

"""
    delta_default(alg) -> Real

Default value of the `delta` option for `alg`.

`delta` is the SDE-specific mixing parameter of the adaptive error estimate: it
weighs the drift (deterministic) error contribution against the diffusion
(stochastic) one when the two are combined into a single `EEst`. A value of `1`
uses the drift estimate at full strength.

Solver subpackages override this for algorithms whose published adaptivity scheme
prescribes a different weighting.
"""
delta_default(alg) = 1 // 1

function OrdinaryDiffEqCore.gamma_default(
        alg::Union{StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm}
    )
    return isadaptive(alg) ? 9 // 10 : 0
end

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
"""
    get_current_alg_order(alg, cache) -> Real

Strong order of the algorithm that is currently being stepped with.

For a plain algorithm this is just `alg_order(alg)`. For a
`StochasticCompositeAlgorithm` it is the order of the member algorithm selected by
`cache.current`, so that adaptivity uses the order of the method that actually took
the step rather than the order of the composite as a whole.
"""
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

"""
    alg_compatible(prob, alg) -> Bool

Whether `alg` can solve `prob`.

`solve` consults this trait before setting up an integrator and errors with a
descriptive message when it returns `false`. Solver subpackages add methods for the
problem classes their algorithm supports, for example a Milstein-type method that
requires diagonal noise:

```julia
alg_compatible(prob::SciMLBase.AbstractSDEProblem, alg::RKMil) = is_diagonal_noise(prob)
```

For a `StochasticCompositeAlgorithm` the result is the maximum over its members, and
for a `JumpProblem` compatibility additionally requires that the algorithm
[`supports_regular_jumps`](@ref) whenever the problem carries a `RegularJump`.
"""
function alg_compatible(
        prob, alg::Union{
            StochasticDiffEqAlgorithm, StochasticDiffEqRODEAlgorithm,
        }
    )
    return true
end
alg_compatible(prob, alg::StochasticDiffEqAlgorithm) = false

function alg_compatible(
        prob::SciMLBase.AbstractSDEProblem,
        alg::Union{
            StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm,
        }
    )
    return max((alg_compatible(prob, a) for a in alg.algs)...)
end

"""
    supports_regular_jumps(alg) -> Bool

Whether `alg` can integrate a `JumpProblem` that carries a `RegularJump`.

Regular jumps are leaped over with a Poisson count per step rather than simulated
one event at a time, which only some SDE integrators implement. The default is
`false`; `EM` and `ImplicitEM` override it to `true` in their subpackages.
[`alg_compatible`](@ref) uses this trait to reject a `RegularJump` problem paired
with an algorithm that cannot handle it.
"""
supports_regular_jumps(alg) = false

# JumpProblem compatibility defaults
function alg_compatible(prob::JumpProblem, alg::StochasticDiffEqAlgorithm)
    return alg_compatible(prob.prob, alg) &&
        (supports_regular_jumps(alg) || prob.regular_jump === nothing) &&
        prob.prob isa SciMLBase.AbstractSDEProblem
end

function alg_compatible(
        prob::JumpProblem,
        alg::Union{StochasticDiffEqJumpAdaptiveAlgorithm, StochasticDiffEqJumpAlgorithm}
    )
    return prob.prob isa DiscreteProblem
end

"""
    alg_needs_extra_process(alg) -> Bool

Whether `alg` needs a second noise process `ΔZ` alongside the Brownian increment `ΔW`.

Higher order schemes (for example the Rößler SRA/SRI families) require an auxiliary
independent process to approximate the extra stochastic integrals appearing in their
order conditions. When this trait is `true` the integrator allocates and resizes the
`Z` process in addition to `W`, so solver subpackages must set it for any algorithm
that reads `integrator.W.dZ`.

For a `StochasticCompositeAlgorithm` the result is the maximum over its members, so
the extra process is present if any member needs it.
"""
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

function OrdinaryDiffEqCore.get_current_alg_autodiff(alg::StochasticDiffEqCompositeAlgorithm, cache)
    return OrdinaryDiffEqCore.alg_autodiff(alg.algs[cache.current])
end

function OrdinaryDiffEqCore.get_chunksize(
        alg::Union{
            StochasticDiffEqNewtonAlgorithm,
            StochasticDiffEqNewtonAdaptiveAlgorithm,
            StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
            StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
        }
    )
    return OrdinaryDiffEqCore._get_fwd_chunksize(typeof(alg.autodiff))
end

OrdinaryDiffEqCore.standardtag(
    alg::Union{
        StochasticDiffEqNewtonAdaptiveAlgorithm,
        StochasticDiffEqNewtonAlgorithm,
    }
) = true

OrdinaryDiffEqCore.alg_difftype(
    alg::Union{
        StochasticDiffEqNewtonAdaptiveAlgorithm,
        StochasticDiffEqNewtonAlgorithm,
    }
) = OrdinaryDiffEqCore._get_fdtype(alg.autodiff)

OrdinaryDiffEqCore.concrete_jac(
    alg::Union{
        StochasticDiffEqNewtonAdaptiveAlgorithm,
        StochasticDiffEqNewtonAlgorithm,
    }
) = alg.concrete_jac

"""
    alg_mass_matrix_compatible(alg) -> Bool

Whether `alg` can solve a problem with a non-identity mass matrix.

`solve` errors when a mass matrix is supplied to an algorithm for which this is
`false`. The implicit theta-method solvers accept a mass matrix only in the
configurations for which the discretization stays consistent (symplectic, or
`theta == 1`) and throw a descriptive error otherwise.
"""
alg_mass_matrix_compatible(alg::StochasticDiffEqAlgorithm) = false
alg_mass_matrix_compatible(alg::StochasticDiffEqRODEAlgorithm) = false

"""
    alg_can_repeat_jac(alg) -> Bool

Whether `alg` may reuse a Jacobian across a rejected-and-repeated step.

When `true` the integrator keeps the factorized `W` matrix after a step rejection
instead of recomputing it, which is the common case. Algorithms whose stage
structure changes between attempts must override this to `false`.
"""
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

"""
    is_split_step(alg) -> Bool

Whether `alg` is a split-step method, i.e. one that advances the drift to an
intermediate state before applying the diffusion rather than adding both
contributions to the same base point.

The split-step solvers (`ISSEM`, `ISSEulerHeun`) override this to `true`; the
integrator uses it to select the matching residual and error-estimate paths.
"""
is_split_step(::StochasticDiffEqAlgorithm) = false

"""
    alg_stability_size(alg) -> Real

Radius of the (deterministic) stability region of `alg` along the negative real axis.

This is the scale used by the automatic stiffness detection of
[`AutoSwitch`](@ref): the stiffness ratio is `abs(eigen_est * dt / alg_stability_size(alg))`,
so an algorithm participating in a stiffness-switching composite must report a
nonzero value. The default of `0` means "not characterized", and is overridden per
algorithm in the solver subpackages.
"""
alg_stability_size(alg::StochasticDiffEqAlgorithm) = 0 # default, overridden per-alg

# is_composite_algorithm trait is defined in OrdinaryDiffEqCore and extended in
# integrator_utils.jl for SDE composite algorithm types.
"""
    unwrap_alg(integrator, is_nlsolve) -> alg

The concrete algorithm that `integrator` should use for the current step.

For a plain algorithm this is `integrator.alg` itself. For a
`StochasticCompositeAlgorithm` it is the currently selected member: the stiff member
when `is_nlsolve` is `true` under [`AutoSwitch`](@ref)-style stiffness switching, and
the member indicated by `integrator.cache.current` otherwise.

`perform_step!` implementations call this instead of reading `integrator.alg`
directly so that they work unchanged inside a composite algorithm.
"""
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

"""
    alg_control_rate(alg) -> Bool

Whether the adaptive step-size controller of `alg` acts on the jump/leap rate rather
than on the usual solution error estimate.

Tau-leaping methods choose `dt` from how much the propensities are allowed to change
over a step, so they override this to `true` and the integrator routes step-size
control through the leap-specific controller instead of the standard `EEst` path.
"""
alg_control_rate(::StochasticDiffEqAlgorithm) = false
alg_control_rate(::StochasticDiffEqRODEAlgorithm) = false
