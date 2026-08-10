module DiffEqBaseMooncakeExt

using DiffEqBase, Mooncake
using DiffEqBase: SciMLBase
using SciMLBase: ADOriginator, ChainRulesOriginator, MooncakeOriginator
import Mooncake: rrule!!, frule!!, CoDual, Dual, zero_fcodual, @is_primitive,
    @from_rrule, @zero_adjoint, @mooncake_overlay, MinimalCtx,
    NoPullback, NoRData, primal, tangent, fdata,
    zero_tangent, instantiate, lazy_zero_rdata

@from_rrule(
    MinimalCtx,
    Tuple{
        typeof(DiffEqBase.solve_up),
        DiffEqBase.AbstractDEProblem,
        Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm},
        Any,
        Any,
        Any,
    },
    true,
)

# Dispatch for auto-alg
@from_rrule(
    MinimalCtx,
    Tuple{
        typeof(DiffEqBase.solve_up),
        DiffEqBase.AbstractDEProblem,
        Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm},
        Any,
        Any,
    },
    true,
)

# promote_f can wrap f in a FunctionWrappersWrapper for performance. Constructing that
# wrapper hits raw :cfunction IR that Mooncake's forward-mode interpreter can't handle
# (SciML/SciMLSensitivity.jl#1427), so this treats the whole call as an opaque primitive
# instead: f/p's tangent is forwarded when returned unchanged, zeroed when f gets wrapped.
# The rare AutoDePSpecialize opaque-p branch isn't handled and errors loudly instead of
# silently returning a wrong gradient.
#
# Only covers the Val{false} method (algorithms that don't need ForwardDiff internally,
# e.g. Tsit5/Verner). The Val{true} method used by algorithms like Rosenbrock/BDF wraps
# jac/tgrad/f separately and is not covered here.
@is_primitive MinimalCtx Tuple{
    typeof(DiffEqBase.promote_f), Any, Val, Any, Any, Any, Val{false}, Val,
}

function _unsupported_promote_f_opaque_p()
    return error(
        "Mooncake differentiation through DiffEqBase.promote_f's AutoDePSpecialize/" *
            "opaque-p branch is not yet supported: `p` is repacked via " *
            "RespecializeParams.pack_auto, which has no derivative rule here.",
    )
end

# `f_out === f_primal` is not just "did promote_f wrap f": promote_f also rebuilds f via
# `@set f.jac_prototype = similar(...)` whenever jac_prototype is set, regardless of
# whether wrapping happens, so f_out can be a fresh (non-identical) but still-unwrapped
# object. Either way, zeroing f_out's tangent below is only exact if f_primal itself had
# nothing to differentiate; otherwise fail loud rather than silently drop a real gradient.
function _check_promote_f_input_has_no_tangent(f_primal)
    Mooncake.tangent_type(typeof(f_primal)) === Mooncake.NoTangent && return nothing
    return error(
        "Mooncake differentiation through DiffEqBase.promote_f is not supported when " *
            "the input `f` (tangent_type = $(Mooncake.tangent_type(typeof(f_primal)))) " *
            "carries its own differentiable state and promote_f returns a new object " *
            "for it (e.g. via wrapping, or via the jac_prototype eltype-promotion path).",
    )
end

function rrule!!(
        pf::CoDual{typeof(DiffEqBase.promote_f)},
        f::CoDual, specialize::CoDual{<:Val}, u0::CoDual, p::CoDual, t::CoDual,
        wrapdiff::CoDual{Val{false}}, cs::CoDual{<:Val},
    )
    f_primal, p_primal = primal(f), primal(p)
    f_out, p_out = DiffEqBase.promote_f(
        f_primal, primal(specialize), primal(u0), p_primal, primal(t), Val(false),
        primal(cs),
    )
    p_out === p_primal || _unsupported_promote_f_opaque_p()
    f_out_is_identity = f_out === f_primal
    f_out_is_identity || _check_promote_f_input_has_no_tangent(f_primal)
    f_out_fdata = f_out_is_identity ? f.dx : fdata(zero_tangent(f_out))
    y = CoDual((f_out, p_out), (f_out_fdata, p.dx))

    # Zero rdata for the argument slots that never reach the outputs (callee, the Val
    # args, and u0/t, whose types but not values matter here). Uses lazy_zero_rdata
    # instead of a hardcoded NoRData() since e.g. t::Float64 has real rdata of its own.
    lazy_pf, lazy_specialize, lazy_u0, lazy_t, lazy_wrapdiff, lazy_cs = map(
        lazy_zero_rdata,
        (primal(pf), primal(specialize), primal(u0), primal(t), primal(wrapdiff), primal(cs)),
    )
    lazy_f = f_out_is_identity ? nothing : lazy_zero_rdata(f_primal)

    # dy's shape depends on whether f_out/p_out carry real rdata: mutable outputs (e.g.
    # array p_out) get their gradient via fdata mutation instead, so dy collapses from a
    # 2-tuple down to a single rdata, or a bare NoRData() if both sides are trivial.
    RDf = Mooncake.rdata_type(Mooncake.tangent_type(typeof(f_out)))
    RDp = Mooncake.rdata_type(Mooncake.tangent_type(typeof(p_out)))
    function promote_f_pb!!(dy)
        df_out, dp_out = if RDf === NoRData && RDp === NoRData
            NoRData(), NoRData()
        elseif RDf === NoRData
            NoRData(), dy
        elseif RDp === NoRData
            dy, NoRData()
        else
            dy
        end
        df = f_out_is_identity ? df_out : instantiate(lazy_f)
        return (
            instantiate(lazy_pf), df, instantiate(lazy_specialize), instantiate(lazy_u0),
            dp_out, instantiate(lazy_t), instantiate(lazy_wrapdiff), instantiate(lazy_cs),
        )
    end
    return y, promote_f_pb!!
end

function frule!!(
        ::Dual{typeof(DiffEqBase.promote_f)},
        f::Dual, specialize::Dual{<:Val}, u0::Dual, p::Dual, t::Dual,
        ::Dual{Val{false}}, cs::Dual{<:Val},
    )
    f_primal, p_primal = primal(f), primal(p)
    f_out, p_out = DiffEqBase.promote_f(
        f_primal, primal(specialize), primal(u0), p_primal, primal(t), Val(false),
        primal(cs),
    )
    p_out === p_primal || _unsupported_promote_f_opaque_p()
    f_out_is_identity = f_out === f_primal
    f_out_is_identity || _check_promote_f_input_has_no_tangent(f_primal)
    f_out_tangent = f_out_is_identity ? tangent(f) : zero_tangent(f_out)
    return Dual((f_out, p_out), (f_out_tangent, tangent(p)))
end

end
