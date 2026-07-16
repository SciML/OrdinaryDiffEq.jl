## SciMLBase Trait Definitions
function SciMLBase.isautodifferentiable(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    return true
end
function SciMLBase.allows_arbitrary_number_types(
        alg::Union{
            OrdinaryDiffEqAlgorithm, DAEAlgorithm,
        }
    )
    return true
end
function SciMLBase.allowscomplex(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    return true
end
# Check if an algorithm's autodiff setting indicates ForwardDiff usage.
# This avoids calling alg_autodiff (defined in OrdinaryDiffEqDifferentiation)
# and handles algorithms without an autodiff field (e.g. ETD2, SplitEuler).
function _autodiff_is_forward(alg)
    hasfield(typeof(alg), :autodiff) || return false
    ad = alg.autodiff
    ad == Val(true) && return true
    ad isa AutoForwardDiff && return true
    ad isa AutoSparse && return dense_ad(ad) isa AutoForwardDiff
    return false
end

function SciMLBase.forwarddiffs_model(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm,
            DAEAlgorithm,
            OrdinaryDiffEqImplicitAlgorithm, ExponentialAlgorithm,
        }
    )
    return _autodiff_is_forward(alg)
end
function SciMLBase.forwarddiffs_model(alg::CompositeAlgorithm)
    return any(_autodiff_is_forward, alg.algs)
end

SciMLBase.forwarddiffs_model_time(alg::RosenbrockAlgorithm) = true

function SciMLBase.forwarddiff_chunksize(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm,
            OrdinaryDiffEqImplicitAlgorithm,
            DAEAlgorithm,
            OrdinaryDiffEqExponentialAlgorithm,
            OrdinaryDiffEqAdaptiveExponentialAlgorithm,
        }
    )
    hasfield(typeof(alg), :autodiff) || return Val{0}()
    return _get_fwd_chunksize(typeof(alg.autodiff))
end

SciMLBase.allows_late_binding_tstops(::OrdinaryDiffEqAlgorithm) = true
SciMLBase.allows_late_binding_tstops(::DAEAlgorithm) = true

SciMLBase.supports_solve_rng(
    ::SciMLBase.AbstractODEProblem,
    ::OrdinaryDiffEqAlgorithm,
) = true

SciMLBase.supports_solve_rng(
    ::SciMLBase.AbstractDAEProblem,
    ::DAEAlgorithm,
) = true

# isadaptive is defined below.

## OrdinaryDiffEq Internal Traits

"""
    isfsal(alg) -> Bool

Return whether `alg` is a FSAL ("first same as last") method, i.e. the last stage
derivative of one step equals the first of the next so it can be reused. Explicit
RK methods are FSAL by default (`true`); the integrator reuses `fsallast` as the
next `fsalfirst` accordingly.
"""
isfsal(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = true
isfsal(tab::DiffEqBase.ExplicitRKTableau) = tab.fsal

# isfsal(alg::CompositeAlgorithm) = isfsal(alg.algs[alg.current])
# Pseudo Non-FSAL
#isfsal(alg::RKM) = false

"""
    isfirk(alg) -> Bool

Return whether `alg` is a fully-implicit Runge–Kutta (FIRK) method
(`false` by default).
"""
isfirk(alg) = false

get_current_isfsal(alg, cache) = isfsal(alg)

"""
    dt_required(alg) -> Bool

Return whether `alg` requires the user to supply a `dt` (`true` by default). Only
fully-adaptive or callback-driven algorithms may relax this.
"""
dt_required(alg) = true

"""
    isdiscretealg(alg) -> Bool

Return whether `alg` is a discrete-time / map-iteration algorithm rather than a
continuous ODE solver (`false` by default).
"""
isdiscretealg(alg) = false

"""
    alg_stability_size(alg) -> Real

Return the (real-axis) linear stability region size of `alg`, used by
stabilized/auto-switching heuristics to decide whether an explicit method is
stable for the current step. Solver sublibraries define it per algorithm.
"""
function alg_stability_size end
"""
    has_stiff_interpolation(alg) -> Bool

Return whether `alg` provides a special interpolant for its stiff branch
(`false` by default).
"""
has_stiff_interpolation(alg) = false

# evaluates f(t[i])
_eval_index(f::F, t::Tuple{A}, _) where {F, A} = f(t[1])
function _eval_index(f::F, t::Tuple{A, Vararg}, i) where {F, A}
    return if i == 1
        f(t[1])
    else
        _eval_index(f, Base.tail(t), i - 1)
    end
end

function get_current_isfsal(alg::CompositeAlgorithm, cache)
    return _eval_index(isfsal, alg.algs, cache.current)::Bool
end

all_fsal(alg, cache) = isfsal(alg)
all_fsal(alg::CompositeAlgorithm, cache) = _all_fsal(alg.algs)

@generated function _all_fsal(algs::T) where {T <: Tuple}
    ex = Expr(
        :tuple, map(1:length(T.types)) do i
            :(isfsal(algs[$i]))
        end...
    )
    return :(all($ex))
end

"""
    issplit(alg) -> Bool

Return whether `alg` treats the RHS as a split function (e.g. IMEX `f = f1 + f2`).
`false` by default.
"""
issplit(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
issplit(alg::StochasticDiffEqAlgorithm) = false

function _composite_beta1_default(
        algs::Tuple{T1, T2}, current, ::Val{QT},
        beta2
    ) where {T1, T2, QT}
    if current == 1
        return QT(beta1_default(algs[1], beta2))
    else
        return QT(beta1_default(algs[2], beta2))
    end
end

@generated function _composite_beta1_default(
        algs::T, current, ::Val{QT},
        beta2
    ) where {T <: Tuple, QT}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(
            expr.args, quote
                if current == $i
                    return QT(beta1_default(algs[$i], beta2))
                end
            end
        )
    end
    return expr
end

function _composite_beta2_default(
        algs::Tuple{T1, T2}, current,
        ::Val{QT}
    ) where {T1, T2, QT}
    if current == 1
        return QT(beta2_default(algs[1]))
    else
        return QT(beta2_default(algs[2]))
    end
end

@generated function _composite_beta2_default(
        algs::T, current,
        ::Val{QT}
    ) where {T <: Tuple, QT}
    expr = Expr(:block)
    for i in 1:length(T.types)
        push!(
            expr.args, quote
                if current == $i
                    return QT(beta2_default(algs[$i]))
                end
            end
        )
    end
    return expr
end

"""
    fsal_typeof(alg, rate_prototype)

Return the type used to store the FSAL derivative for `alg` given a
`rate_prototype`. Defaults to `typeof(rate_prototype)`; overridden by algorithms
whose FSAL slot differs from the ordinary rate type.
"""
function fsal_typeof(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, rate_prototype)
    return typeof(rate_prototype)
end

function fsal_typeof(alg::CompositeAlgorithm, rate_prototype)
    fsal = map(x -> fsal_typeof(x, rate_prototype), alg.algs)
    @assert length(unique(fsal)) == 1 "`fsal_typeof` must be consistent"
    return fsal[1]
end

"""
    isimplicit(alg) -> Bool

Return whether `alg` solves an implicit (nonlinear) system each step. `false` for
explicit methods; `true` for implicit ones. For a [`CompositeAlgorithm`](@ref) it
is `true` if any constituent is implicit.
"""
isimplicit(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isimplicit(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
isimplicit(alg::OrdinaryDiffEqImplicitAlgorithm) = true
isimplicit(alg::CompositeAlgorithm) = any(isimplicit.(alg.algs))

"""
    isdtchangeable(alg) -> Bool

Return whether `alg` allows the step size `dt` to change between steps (`true` by
default). Multistep-style methods that require a fixed grid set this to `false`.
"""
isdtchangeable(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = true
isdtchangeable(alg::CompositeAlgorithm) = all(isdtchangeable.(alg.algs))
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
isdtchangeable(alg) = true

"""
    ismultistep(alg) -> Bool

Return whether `alg` is a multistep method (uses solution history from more than
the previous step). `false` by default.
"""
ismultistep(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
ismultistep(alg::CompositeAlgorithm) = any(ismultistep.(alg.algs))

"""
    isadaptive(alg) -> Bool

Return whether `alg` performs adaptive step-size control. `true` for algorithms
subtyping an `…Adaptive…` abstract type, `false` otherwise.
"""
isadaptive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
isadaptive(alg::OrdinaryDiffEqAdaptiveAlgorithm) = true
isadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = all(isadaptive.(alg.algs))
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
isadaptive(alg) = false

"""
    has_special_newton_error(alg) -> Bool

Return whether `alg` supplies its own Newton-iteration error estimate rather than
the generic one (`false` by default).
"""
has_special_newton_error(alg) = false

anyadaptive(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = isadaptive(alg)
anyadaptive(alg::OrdinaryDiffEqCompositeAlgorithm) = any(isadaptive, alg.algs)
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
anyadaptive(alg) = isadaptive(alg)

"""
    has_dtnew_modification(alg) -> Bool

Return whether `alg` post-processes the controller's proposed `dtnew` via a
`dtnew_modification` hook (`false` by default).
"""
has_dtnew_modification(alg) = false
dtnew_modification(integrator, alg, dtnew) = dtnew

# Whether the algorithm's alg_cache / perform_step! handle `u === nothing` directly.
# Default: false, and __init coerces a null u0 to Float64[] so generic cache code
# (zero(u), similar(u), broadcasts) still works. Solvers that want to preserve the
# `nothing` signal — e.g., to bypass nonlinear solvers or skip allocation — override
# to true.
"""
    allows_null_u0(alg) -> Bool

Return whether `alg` supports an empty (zero-length) initial condition `u0`
(`false` by default).
"""
allows_null_u0(alg) = false

# Whether an algorithm uses a posteriori dt estimates (always accepts, then picks next dt).
# Default is false. CaoTauLeaping overrides to true.
"""
    isaposteriori(alg) -> Bool

Return whether `alg` uses a-posteriori (rather than embedded) error estimation
(`false` by default).
"""
isaposteriori(alg) = false

"""
    isautoswitch(alg) -> Bool

Return whether `alg` is a composite algorithm whose choice function is an
[`AutoSwitch`](@ref) stiffness detector.
"""
isautoswitch(alg) = false
isautoswitch(alg::CompositeAlgorithm) = alg.choice_function isa AutoSwitch

"""
    only_diagonal_mass_matrix(alg) -> Bool

Return whether `alg` only supports diagonal (not general) mass matrices
(`false` by default). Used to decide FSAL reevaluation.
"""
only_diagonal_mass_matrix(alg) = false
"""
    isdp8(alg) -> Bool

Return whether `alg` is the `DP8` method. Used to special-case its FSAL / dense
handling (`false` by default).
"""
isdp8(alg) = false
"""
    isdefaultalg(alg) -> Bool

Return whether `alg` is the automatic default algorithm wrapper. Used to route the
stiffness auto-switching logic through its specialized default path.
"""
isdefaultalg(alg) = false

"""
    qmin_default(alg) -> Real

Return the algorithm-specific lower bound for the per-step shrink factor.

`resolve_basic` uses this developer extension point to fill
[`CommonControllerOptions`](@ref) when the user did not pass `qmin`.

# Arguments

- `alg`: The solver algorithm whose controller defaults are being resolved.

# Returns

- `Real`: The default lower bound. Adaptive algorithms use `1 // 5`; non-adaptive
    algorithms use `0`.
"""
function qmin_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    return isadaptive(alg) ? 1 // 5 : 0
end
qmin_default(alg::CompositeAlgorithm) = maximum(qmin_default.(alg.algs))
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
qmin_default(alg) = 1 // 5

"""
    qmax_default(alg) -> Real

Return the algorithm-specific upper bound for the per-step growth factor.

`resolve_basic` uses this developer extension point to fill
[`CommonControllerOptions`](@ref) when the user did not pass `qmax`.

# Arguments

- `alg`: The solver algorithm whose controller defaults are being resolved.

# Returns

- `Real`: The default upper bound. The generic default is `10`.
"""
qmax_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 10
qmax_default(alg::CompositeAlgorithm) = minimum(qmax_default.(alg.algs))
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
qmax_default(alg) = 10

"""
    get_chunksize(alg) -> Val

Return the ForwardDiff chunk size configured on `alg`'s autodiff choice, as a
`Val`. `Val(0)` means "let ForwardDiff choose".
"""
function get_chunksize(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have a chunk size defined.")
end

"""
    _get_fwd_chunksize(AD) -> Val

Return, as a `Val`, the ForwardDiff chunk size encoded in the `AutoForwardDiff`
type/instance `AD` (`Val(0)` when unspecified).
"""
_get_fwd_chunksize(::Type{<:AutoForwardDiff{nothing}}) = Val(0)
_get_fwd_chunksize(::Type{<:AutoForwardDiff{CS}}) where {CS} = Val(CS)
"""
    _get_fwd_chunksize_int(AD) -> Int

Return the ForwardDiff chunk size encoded in the `AutoForwardDiff` type/instance
`AD` as a plain `Int` (`0` when unspecified).
"""
_get_fwd_chunksize_int(::Type{<:AutoForwardDiff{nothing}}) = 0
_get_fwd_chunksize_int(::Type{<:AutoForwardDiff{CS}}) where {CS} = CS
_get_fwd_chunksize(AD) = Val(0)
_get_fwd_chunksize_int(AD) = 0
_get_fwd_chunksize_int(::AutoForwardDiff{nothing}) = 0
_get_fwd_chunksize_int(::AutoForwardDiff{CS}) where {CS} = CS
_get_fwd_tag(::AutoForwardDiff{CS, T}) where {CS, T} = T

"""
    _get_fdtype(AD)

Return the finite-difference type parameter (e.g. `Val{:forward}`) of an
`AutoFiniteDiff` type/instance `AD`.
"""
_get_fdtype(::AutoFiniteDiff{T1}) where {T1} = T1
_get_fdtype(::Type{<:AutoFiniteDiff{T1}}) where {T1} = T1

function get_chunksize(
        alg::Union{
            OrdinaryDiffEqExponentialAlgorithm,
            OrdinaryDiffEqAdaptiveExponentialAlgorithm,
            OrdinaryDiffEqImplicitAlgorithm,
            OrdinaryDiffEqAdaptiveImplicitAlgorithm,
            DAEAlgorithm,
        }
    )
    hasfield(typeof(alg), :autodiff) || return Val{0}()
    return _get_fwd_chunksize(typeof(alg.autodiff))
end

function get_chunksize_int(alg::OrdinaryDiffEqAlgorithm)
    error("This algorithm does not have a chunk size defined.")
end

function get_chunksize_int(
        alg::Union{
            OrdinaryDiffEqExponentialAlgorithm,
            OrdinaryDiffEqAdaptiveExponentialAlgorithm,
            OrdinaryDiffEqImplicitAlgorithm,
            OrdinaryDiffEqAdaptiveImplicitAlgorithm,
            DAEAlgorithm,
        }
    )
    hasfield(typeof(alg), :autodiff) || return 0
    return _get_fwd_chunksize_int(typeof(alg.autodiff))
end

# get_chunksize(alg::CompositeAlgorithm) = get_chunksize(alg.algs[alg.current_alg])

"""
    alg_autodiff(alg)

Return the automatic-differentiation choice (an ADTypes.jl `AbstractADType`) that
`alg` uses to build Jacobians. Implicit-solver sublibraries define this for their
algorithms; see also [`get_current_alg_autodiff`](@ref) for composite algorithms.
"""
function alg_autodiff end

# Linear Exponential doesn't have any of the AD stuff
function DiffEqBase.prepare_alg(
        alg::OrdinaryDiffEqLinearExponentialAlgorithm,
        u0::AbstractArray,
        p, prob
    )
    return alg
end

function DiffEqBase.prepare_alg(alg::CompositeAlgorithm, u0, p, prob)
    algs = map(a -> DiffEqBase.prepare_alg(a, u0, p, prob), alg.algs)
    cf = alg.choice_function
    if cf isa AutoSwitch
        nonstiffalg = _prepare_autoswitch_alg(cf.nonstiffalg, u0, p, prob)
        stiffalg = _prepare_autoswitch_alg(cf.stiffalg, u0, p, prob)
        cf = AutoSwitch(
            nonstiffalg, stiffalg,
            cf.maxstiffstep, cf.maxnonstiffstep,
            cf.nonstifftol, cf.stifftol,
            cf.dtfac, cf.stiffalgfirst, cf.switch_max
        )
    end
    return CompositeAlgorithm(algs, cf)
end

_prepare_autoswitch_alg(alg, u0, p, prob) = DiffEqBase.prepare_alg(alg, u0, p, prob)
function _prepare_autoswitch_alg(algs::Tuple, u0, p, prob)
    return map(a -> DiffEqBase.prepare_alg(a, u0, p, prob), algs)
end

"""
    has_autodiff(alg) -> Bool

Return whether `alg` carries an `autodiff` field configuring Jacobian
differentiation (`false` for explicit methods).
"""
has_autodiff(alg::OrdinaryDiffEqAlgorithm) = false
function has_autodiff(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm, OrdinaryDiffEqImplicitAlgorithm,
            CompositeAlgorithm, DAEAlgorithm,
        }
    )
    return true
end
# ExponentialAlgorithm subtypes may or may not have an autodiff field
# (e.g. ETD2, SplitEuler, LinearExponential do not)
has_autodiff(alg::ExponentialAlgorithm) = hasfield(typeof(alg), :autodiff)

# end

# alg_autodiff(alg::CompositeAlgorithm) = alg_autodiff(alg.algs[alg.current_alg])
"""
    get_current_alg_autodiff(alg, cache)

Return the autodiff choice active for the current step. Equals
[`alg_autodiff`](@ref)`(alg)` for simple algorithms; for a
[`CompositeAlgorithm`](@ref) it selects the currently-active constituent via
`cache.current`.
"""
get_current_alg_autodiff(alg, cache) = alg_autodiff(alg)
function get_current_alg_autodiff(alg::CompositeAlgorithm, cache)
    return _eval_index(alg_autodiff, alg.algs, cache.current)::Bool
end

"""
    alg_difftype(alg)

Return the finite-difference type (e.g. `Val{:forward}`) configured on `alg`'s
`AutoFiniteDiff` autodiff choice, used when differentiating by finite differences.
"""
function alg_difftype(
        alg::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm,
            OrdinaryDiffEqImplicitAlgorithm,
            OrdinaryDiffEqExponentialAlgorithm,
            OrdinaryDiffEqAdaptiveExponentialAlgorithm,
            DAEAlgorithm,
        }
    )
    hasfield(typeof(alg), :autodiff) || return Val{:forward}
    return _get_fdtype(alg.autodiff)
end

"""
    standardtag(alg) -> Bool

Return whether `alg` uses the standard ForwardDiff tagging (a custom tag type for
Dual numbers) when building Jacobians. Used by the differentiation machinery to
pick the AD config.
"""
standardtag(
    alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm,
        OrdinaryDiffEqImplicitAlgorithm,
        OrdinaryDiffEqExponentialAlgorithm,
        OrdinaryDiffEqAdaptiveExponentialAlgorithm,
        DAEAlgorithm,
    }
) = true

"""
    concrete_jac(alg)

Return the `concrete_jac` setting of `alg`: `true`/`false` forces whether a
concrete Jacobian matrix is materialized, `nothing` lets the solver decide (e.g.
based on the chosen linear solver).
"""
concrete_jac(
    alg::Union{
        OrdinaryDiffEqAdaptiveImplicitAlgorithm,
        OrdinaryDiffEqImplicitAlgorithm,
        OrdinaryDiffEqExponentialAlgorithm,
        OrdinaryDiffEqAdaptiveExponentialAlgorithm,
        DAEAlgorithm,
    }
) = alg.concrete_jac

"""
    alg_extrapolates(alg) -> Bool

Return whether `alg` needs an extrapolated initial guess for its implicit/predictor
stage (`false` by default). Algorithms that do set this to `true` so the
integrator computes `uprev2`/extrapolant state for them.
"""
alg_extrapolates(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
alg_extrapolates(alg::CompositeAlgorithm) = any(alg_extrapolates.(alg.algs))
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
alg_extrapolates(alg) = false
"""
    alg_order(alg) -> Int

Return the order of accuracy of `alg`. Solver sublibraries define this for each
concrete algorithm; for a [`CompositeAlgorithm`](@ref) it is the maximum over the
constituent algorithms.
"""
function alg_order(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    error("Order is not defined for this algorithm")
end
alg_order(alg::CompositeAlgorithm) = maximum(alg_order, alg.algs)

"""
    get_current_alg_order(alg, cache) -> Int

Return the order currently in effect for `alg` given its `cache`. Equal to
[`alg_order`](@ref) for fixed-order methods; for variable-order methods
(Adams/BDF) it reads the order stored on the cache.
"""
function get_current_alg_order(
        alg::Union{
            OrdinaryDiffEqAlgorithm, DAEAlgorithm,
            StochasticDiffEqAlgorithm,
            StochasticDiffEqRODEAlgorithm,
        }, cache
    )
    return alg_order(alg)
end
function get_current_alg_order(alg::CompositeAlgorithm, cache)
    return _eval_index(alg_order, alg.algs, cache.current)::Int
end

get_current_alg_order(alg::OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm, cache) = cache.order
"""
    get_current_adaptive_order(alg, cache) -> Int

Return the order used for the current step's adaptive error estimate given `alg`
and its `cache` (order-aware analogue of [`alg_adaptive_order`](@ref) for
variable-order methods).
"""
function get_current_adaptive_order(alg::OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm, cache)
    return cache.order
end

#alg_adaptive_order(alg::OrdinaryDiffEqAdaptiveAlgorithm) = error("Algorithm is adaptive with no order")
function get_current_adaptive_order(alg, cache)
    return alg_adaptive_order(alg)
end
function get_current_adaptive_order(alg::CompositeAlgorithm, cache)
    return _eval_index(alg_adaptive_order, alg.algs, cache.current)::Int
end

"""
    alg_maximum_order(alg) -> Int

Return the maximum order the algorithm can attain. Equal to [`alg_order`](@ref)
for fixed-order methods; for composite algorithms it is the maximum over the
constituents.
"""
alg_maximum_order(alg) = alg_order(alg)
alg_maximum_order(alg::CompositeAlgorithm) = maximum(alg_order(x) for x in alg.algs)

"""
    alg_adaptive_order(alg) -> Int

Return the order used for the adaptive error estimate of `alg`. The generic
fallback is `alg_order(alg) - 1`; it is deliberately conservative because it
tracks the realized error better than the embedded-estimate order.
"""
alg_adaptive_order(alg) = alg_order(alg) - 1

# this is actually incorrect and is purposefully decreased as this tends
# to track the real error much better

"""
    default_controller(::Type{QT}, alg)

Return the step-size controller used by `alg` when the user does not pass
one explicitly via `solve(prob, alg; controller = ...)`. `QT` is the
scalar type used internally for `q`, `dt`, and the step-size factors —
usually `Float64`. The generic fallback is `PIController(QT, alg)`;
override for new algorithm types whose default should differ
(e.g. BDF uses `BDFController`, JVODE uses `JVODEController`).
"""
function default_controller(QT, alg)
    return PIController(QT, alg)
end

function default_controller(QT, alg::OrdinaryDiffEqCompositeAlgorithm)
    return CompositeController(
        map(alg -> default_controller(QT, alg), alg.algs)
    )
end

function _digest_beta1_beta2(alg, cache, ::Val{QT}, _beta1, _beta2) where {QT}
    if alg isa OrdinaryDiffEqCompositeAlgorithm
        beta2 = _beta2 === nothing ?
            _composite_beta2_default(alg.algs, cache.current, Val(QT)) : _beta2
        beta1 = _beta1 === nothing ?
            _composite_beta1_default(alg.algs, cache.current, Val(QT), beta2) : _beta1
    else
        beta2 = _beta2 === nothing ? beta2_default(alg) : _beta2
        beta1 = _beta1 === nothing ? beta1_default(alg, beta2) : _beta1
    end
    return convert(QT, beta1)::QT, convert(QT, beta2)::QT
end

"""
    beta2_default(alg) -> Real

Return the algorithm-specific default integral gain for [`PIController`](@ref).

# Arguments

- `alg`: The solver algorithm whose PI-controller defaults are being resolved.

# Returns

- `Real`: For adaptive algorithms, `2 // (5 * alg_order(alg))`; otherwise `0`.
"""
function beta2_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    return isadaptive(alg) ? 2 // (5alg_order(alg)) : 0
end

"""
    beta1_default(alg, beta2) -> Real

Return the algorithm-specific default proportional gain for
[`PIController`](@ref).

`beta1` is allowed to depend on the resolved `beta2` value, so algorithms that
override both should normally override `beta2_default` first.

# Arguments

- `alg`: The solver algorithm whose PI-controller defaults are being resolved.
- `beta2`: The resolved integral gain.

# Returns

- `Real`: For adaptive algorithms, `7 // (10 * alg_order(alg))`; otherwise `0`.
"""
function beta1_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, beta2)
    return isadaptive(alg) ? 7 // (10alg_order(alg)) : 0
end

"""
    gamma_default(alg)

Algorithm-specific default for the safety factor on the controller's
predicted dt change. The generic fallback is `9 // 10` for adaptive
algorithms and `0` otherwise.
"""
function gamma_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm})
    return isadaptive(alg) ? 9 // 10 : 0
end
gamma_default(alg::CompositeAlgorithm) = maximum(gamma_default, alg.algs)
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
gamma_default(alg) = isadaptive(alg) ? 9 // 10 : 0

# Whether `PredictiveController.stepsize_controller!` should use the
# plain `gamma` factor without the Newton-iter correction (used by FIRK).
"""
    fac_default_gamma(alg) -> Bool

Return whether the [`PredictiveController`](@ref) should apply the plain `gamma`
safety factor without the Newton-iteration correction for `alg` (`false` by
default; `true` for FIRK methods).
"""
fac_default_gamma(alg) = false

"""
    qsteady_min_default(alg) -> Real

Return the lower edge of the algorithm-specific step-size deadband.

If the proposed step-size factor lies inside
`[qsteady_min, qsteady_max]`, the controller holds `dt` constant.

# Arguments

- `alg`: The solver algorithm whose controller deadband is being resolved.

# Returns

- `Real`: The lower deadband edge. The generic default is `1`.
"""
qsteady_min_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 1

"""
    qsteady_max_default(alg) -> Real

Return the upper edge of the algorithm-specific step-size deadband.

If the proposed step-size factor lies inside
`[qsteady_min, qsteady_max]`, the controller holds `dt` constant.
Adaptive implicit algorithms widen this default to `6 // 5` to reduce
Jacobian recomputation; BDF methods may specialize it further.

# Arguments

- `alg`: The solver algorithm whose controller deadband is being resolved.

# Returns

- `Real`: The upper deadband edge. The generic default is `1`.
"""
qsteady_max_default(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = 1
# Generic fallbacks for non-ODE algorithms (SDE, RODE) calling __init
qsteady_min_default(alg) = 1
qsteady_max_default(alg) = 1
qsteady_max_default(alg::OrdinaryDiffEqAdaptiveImplicitAlgorithm) = 6 // 5
# But don't re-use Jacobian if not adaptive: too risky and cannot pull back
qsteady_max_default(alg::OrdinaryDiffEqImplicitAlgorithm) = isadaptive(alg) ? 1 // 1 : 0

"""
    qmax_first_step_default(alg)

Algorithm-specific default for the looser `qmax` applied to the very
first accepted step (mirrors Sundials CVODE — the initial dt from the
automatic step-size selection is approximate, so a much larger growth
is allowed once). Generic fallback is `10000`.
See https://github.com/SciML/DifferentialEquations.jl/issues/299.
"""
qmax_first_step_default(alg) = 10000

"""
    failfactor_default(alg)

Algorithm-specific default for the post-Newton-failure dt shrink factor
used by [`post_newton_controller!`](@ref). Generic fallback is `2`
(halve `dt` on each failed Newton solve).
"""
failfactor_default(alg) = 2
#TODO
#SciMLBase.nlsolve_default(::QNDF, ::Val{κ}) = 1//2

# SSP coefficients

"""
    ssp_coefficient(alg)

Return the SSP coefficient of the ODE algorithm `alg`. If one time step of size
`dt` with `alg` can be written as a convex combination of explicit Euler steps
with step sizes `cᵢ * dt`, the SSP coefficient is the minimal value of `1/cᵢ`.

# Examples

```julia-repl
julia> ssp_coefficient(SSPRK104())
6
```
"""
ssp_coefficient(alg) = error("$alg is not a strong stability preserving method.")

# We shouldn't do this probably.
#ssp_coefficient(alg::ImplicitEuler) = Inf

"""
    alg_can_repeat_jac(alg) -> Bool

Return whether `alg` may reuse a Jacobian/W factorization across steps rather than
recomputing every step. `true` for adaptive Newton algorithms, `false` otherwise.
"""
alg_can_repeat_jac(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
alg_can_repeat_jac(alg::OrdinaryDiffEqNewtonAdaptiveAlgorithm) = true

"""
    unwrap_alg(integrator, is_stiff)
    unwrap_alg(alg, is_stiff)

Return the concrete algorithm actually driving the current step. For a
non-composite algorithm this is `alg` itself; for a [`CompositeAlgorithm`](@ref)
it selects the active constituent based on the current cache index or, for a
two-member auto-switch pair, on the `is_stiff` flag.
"""
function unwrap_alg(alg::SciMLBase.AbstractDEAlgorithm, is_stiff)
    if !is_composite_algorithm(alg)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        if length(alg.algs) > 2
            return alg.algs[alg.choice_function.current]
        end
        if is_stiff === nothing
            throwautoswitch(alg)
        end
        num = is_stiff ? 2 : 1
        if num == 1
            return alg.algs[1]
        else
            return alg.algs[2]
        end
    else
        error("this dispatch does not support this algorithm right now")
    end
end

function unwrap_alg(integrator, is_stiff)
    alg = integrator.alg
    if !is_composite_algorithm(alg)
        return alg
    elseif alg.choice_function isa AutoSwitchCache
        if length(alg.algs) > 2
            alg.algs[alg.choice_function.current]
        else
            if is_stiff === nothing
                throwautoswitch(alg)
            end
            num = is_stiff ? 2 : 1
            if num == 1
                return alg.algs[1]
            else
                return alg.algs[2]
            end
        end
    else
        return _eval_index(identity, alg.algs, integrator.cache.current)
    end
end

function throwautoswitch(alg)
    throw(ArgumentError("one of $(alg.algs) is not compatible with stiffness-based autoswitching"))
end

# Whether `uprev` is used in the algorithm directly.
"""
    uses_uprev(alg, adaptive::Bool) -> Bool

Return whether `alg` uses the previous step value `uprev` directly (as opposed to
only via the FSAL history). `true` by default; used to decide whether the `uprev`
cache slot must be maintained.
"""
uses_uprev(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}, adaptive::Bool) = true
uses_uprev(alg::OrdinaryDiffEqAdaptiveAlgorithm, adaptive::Bool) = true
# Generic fallback for non-ODE algorithms (SDE, RODE) calling __init
uses_uprev(alg, adaptive) = true


"""
    isWmethod(alg) -> Bool

Return whether `alg` is a W-method, i.e. it remains correct with an inexact
(stale) Jacobian in its `W` matrix (`false` by default). Rosenbrock-W methods set
this `true`.
"""
isWmethod(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

"""
    isesdirk(alg) -> Bool

Return whether `alg` is an ESDIRK (explicit-first-stage singly-diagonally-implicit
Runge–Kutta) method (`false` by default).
"""
isesdirk(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false

"""
    is_mass_matrix_alg(alg) -> Bool

Return whether `alg` supports solving mass-matrix ODEs/DAEs `M u' = f`
(`false` by default; `true` for Rosenbrock and appropriate Newton methods).
"""
is_mass_matrix_alg(alg::Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm}) = false
is_mass_matrix_alg(alg::CompositeAlgorithm) = all(is_mass_matrix_alg, alg.algs)
is_mass_matrix_alg(alg::RosenbrockAlgorithm) = true
is_mass_matrix_alg(alg::NewtonAlgorithm) = isfsal(alg) || !isesdirk(alg)

# All algorithms should be shown using their keyword definition, and not as structs
function Base.show(io::IO, ::MIME"text/plain", alg::OrdinaryDiffEqAlgorithm)
    print(io, String(typeof(alg).name.name), "(;")
    for fieldname in fieldnames(typeof(alg))
        print(io, " ", fieldname, " = ", getfield(alg, fieldname), ",")
    end
    return print(io, ")")
end

# Defaults in the current system: currently opt out DAEAlgorithms until complete
"""
    default_linear_interpolation(alg, prob) -> Bool

Return whether the solver should default to (cheaper) linear interpolation instead
of the algorithm's Hermite/dense interpolant for `alg` on `prob`. `true` for DAEs,
discrete problems, and RODE/SDE problems.
"""
default_linear_interpolation(alg, prob) = alg isa DAEAlgorithm || prob isa DiscreteProblem
# RODE/SDE always uses linear interpolation (no dense output)
default_linear_interpolation(prob::SciMLBase.AbstractRODEProblem, alg) = true
