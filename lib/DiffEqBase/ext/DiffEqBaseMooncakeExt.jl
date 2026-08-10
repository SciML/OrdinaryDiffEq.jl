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

# `promote_f`'s "simple" method (used when the algorithm does not need ForwardDiff
# internally, e.g. Tsit5/Verner): either returns `(f, p)` unchanged (no wrapping needed),
# or wraps `f` in a `FunctionWrappersWrapper` for performance (type-erasing dynamic
# dispatch), passing `p` straight through in the common case, or (only in the rare
# `AutoDePSpecialize` "opaque" branch) repacking it via `RespecializeParams.pack_auto`.
# None of this can be traced through directly in forward mode: constructing a
# `FunctionWrappersWrapper` internally hits raw `:cfunction` IR that Mooncake's
# forward-mode interpreter does not support (SciML/SciMLSensitivity.jl#1427). This treats
# the whole call as an opaque primitive with a hand-derived tangent instead: `f`'s own
# tangent/fdata is forwarded when it's returned unchanged, `NoTangent`/`NoFData` when it's
# wrapped (matching `FunctionWrappersWrapper`'s own declared `NoTangent` tangent type), and
# `p`'s tangent/fdata is forwarded unchanged in the identity case. The
# `AutoDePSpecialize`/opaque-p case isn't handled -- it fails loudly with a clear error
# instead of silently returning a wrong gradient.
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
    f_out_fdata = f_out_is_identity ? f.dx : fdata(zero_tangent(f_out))
    y = CoDual((f_out, p_out), (f_out_fdata, p.dx))

    # Zero rdata for every argument slot that never participates in the two real
    # outputs (`f_out`/`p_out`): the callee itself, `specialize`/`wrapdiff`/`cs` (all
    # `Val`s), and `u0`/`t` (only their *types*, not values, matter here -- they're
    # solely used to build the wrapped function's dispatch signature). Computed lazily
    # via `lazy_zero_rdata`/`instantiate` rather than a hardcoded `NoRData()`: several of
    # these (e.g. `t::Float64`) have a genuinely non-trivial rdata type even though their
    # *value* has zero contribution here, and returning the wrong type for those would
    # break accumulation at the call site.
    lazy_pf, lazy_specialize, lazy_u0, lazy_t, lazy_wrapdiff, lazy_cs = map(
        lazy_zero_rdata,
        (primal(pf), primal(specialize), primal(u0), primal(t), primal(wrapdiff), primal(cs)),
    )
    lazy_f = f_out_is_identity ? nothing : lazy_zero_rdata(f_primal)

    # `dy`'s shape depends on how Mooncake collapses `(rdata_type(f_out),
    # rdata_type(p_out))`: mutable types (e.g. array `p_out`) carry their real gradient
    # via `fdata` mutation, not `rdata`, so their rdata contribution is trivial
    # (`NoRData`). When one or both slots are trivial, the incoming `dy` isn't a 2-tuple:
    # it's the single remaining non-trivial rdata, or a bare `NoRData()` if both are.
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
    f_out_tangent = f_out === f_primal ? tangent(f) : zero_tangent(f_out)
    return Dual((f_out, p_out), (f_out_tangent, tangent(p)))
end

end
