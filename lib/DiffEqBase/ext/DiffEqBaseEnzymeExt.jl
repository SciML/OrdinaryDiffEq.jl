module DiffEqBaseEnzymeExt

@static if isempty(VERSION.prerelease)
    using DiffEqBase
    import DiffEqBase: value
    using Enzyme
    import Enzyme: Const, MixedDuplicated
    using ChainRulesCore
    import SciMLStructures


    @inline function copy_or_reuse(config, val, idx)
        if Enzyme.EnzymeRules.overwritten(config)[idx] && ismutable(val)
            return deepcopy(val)
        else
            return val
        end
    end

    @inline function arg_copy(data, i)
        config, args = data
        copy_or_reuse(config, args[i].val, i + 5)
    end

    # Accumulate a tangent `darg` into a shadow `dval`.
    # When `dval` is a SciMLStructure (e.g. MTKParameters), `darg` may be:
    #   - A tunable gradient vector (from SciMLSensitivity's EnzymeOriginator path)
    #   - Another SciMLStructure
    #   - A broadcastable array
    # In all cases, accumulation goes through the SciMLStructures interface.
    #
    # `diff_tunables` mirrors the sensealg field of the same name and means
    # "differentiate only the Tunable portion." When `true` (the default and
    # the value carried by `QuadratureAdjoint`/`GaussAdjoint`/`InterpolatingAdjoint`
    # adjoints unless the user opted out) only the Tunable slice of a structured
    # `darg` is accumulated. When `false`, `SciMLSensitivity.adjointbackpass`
    # returns a structured cotangent whose gradient contribution may live in
    # non-Tunable fields such as `caches` (e.g. SCC sub-problem buffers feeding
    # `explicitfuns!`), so those fields are walked in as well.
    function _accum_tangent!(dval, darg; diff_tunables::Bool = true)
        if SciMLStructures.isscimlstructure(dval) && !(dval isa AbstractArray)
            if SciMLStructures.isscimlstructure(darg)
                shadow_tunables, _, _ = SciMLStructures.canonicalize(
                    SciMLStructures.Tunable(), dval,
                )
                darg_tunables, _, _ = SciMLStructures.canonicalize(
                    SciMLStructures.Tunable(), darg,
                )
                shadow_tunables .+= darg_tunables
                SciMLStructures.replace!(SciMLStructures.Tunable(), dval, shadow_tunables)
                if !diff_tunables
                    for field in fieldnames(typeof(darg))
                        field === :tunable && continue
                        hasfield(typeof(dval), field) || continue
                        _accum_nested!(getfield(dval, field), getfield(darg, field))
                    end
                end
            elseif darg isa AbstractVector
                shadow_tunables, _, _ = SciMLStructures.canonicalize(
                    SciMLStructures.Tunable(), dval,
                )
                if length(darg) == length(shadow_tunables)
                    shadow_tunables .+= darg
                    SciMLStructures.replace!(
                        SciMLStructures.Tunable(), dval, shadow_tunables,
                    )
                else
                    dval .+= darg
                end
            elseif darg isa NamedTuple
                # Full parameter gradient from a structured adjoint backpass
                # (includes caches and other non-tunable components).
                # Accumulate each matching field into the shadow.
                for field in fieldnames(typeof(darg))
                    src = getfield(darg, field)
                    src === nothing && continue
                    if hasfield(typeof(dval), field)
                        dst = getfield(dval, field)
                        _accum_nested!(dst, src)
                    end
                end
            else
                dval .+= darg
            end
        else
            dval .+= darg
        end
        return nothing
    end

    # Recursively accumulate nested containers (tuples of arrays, etc.)
    function _accum_nested!(dst::AbstractArray, src::AbstractArray)
        dst .+= src
        return nothing
    end
    function _accum_nested!(dst::Tuple, src::Tuple)
        for (d, s) in zip(dst, src)
            _accum_nested!(d, s)
        end
        return nothing
    end
    _accum_nested!(::Any, ::Nothing) = nothing
    _accum_nested!(::Nothing, ::Any) = nothing
    _accum_nested!(::Nothing, ::Nothing) = nothing

    # Build the shadow `ODESolution` for the augmented primal. A plain
    # `Enzyme.make_zero(sol)` recursively allocates fresh zero buffers for every
    # mutable field of `sol`, including `sol.prob.p` and `sol.prob.u0`, which are
    # aliased to the active `p` / `u0` shadows the outer caller is differentiating
    # against. Severing that aliasing means a cotangent written into the returned
    # `sol.prob.p` field by a downstream consumer (e.g. anything reading the
    # solution's parameters) goes into a dangling buffer instead of the buffer
    # the outer Enzyme tape is tracking, silently dropping that gradient
    # contribution.
    #
    # Pre-seed the `make_zero` seen-set so `prob.p` and `prob.u0` map to
    # themselves: recursion into those fields short-circuits via `haskey(seen, …)`,
    # preserving aliasing with the outer shadow. The actual derivative-carrying
    # field (`sol.u`) still gets a fresh zero buffer.
    @inline function _make_solution_zero(sol)
        seen = IdDict()
        _preseed_alias!(seen, sol.prob.p)
        _preseed_alias!(seen, sol.prob.u0)
        return Enzyme.make_zero(Core.Typeof(sol), seen, sol, Val(false))
    end

    @inline _preseed_alias!(::IdDict, ::Nothing) = nothing
    @inline function _preseed_alias!(seen::IdDict, v)
        if ismutable(v)
            seen[v] = v
        end
        return nothing
    end

    # Note these following functions are generally not considered user facing from within Enzyme.
    # They enable additional performance/usability here (e.g. inactive kwargs).
    # Contact wsmoses@ before modifying (and beware their semantics may change without semver).

    Enzyme.EnzymeRules.inactive_kwarg(::typeof(DiffEqBase.solve_up), prob, sensalg::Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}, u0, p, args...; kwargs...) = nothing

    Enzyme.EnzymeRules.has_easy_rule(::typeof(DiffEqBase.solve_up), prob, sensalg::Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}, u0, p, args...; kwargs...) = nothing

    function Enzyme.EnzymeRules.augmented_primal(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, RTA::Union{Type{Duplicated{RT}}, Type{MixedDuplicated{RT}}}, prob,
            sensealg::Union{
                Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm},
            },
            u0, p, args...; kwargs...
        ) where {RT}

        res = DiffEqBase._solve_adjoint(
            copy_or_reuse(config, prob.val, 2), copy_or_reuse(config, sensealg.val, 3),
            copy_or_reuse(config, u0.val, 4), copy_or_reuse(config, p.val, 5),
            SciMLBase.EnzymeOriginator(), ntuple(Base.Fix1(arg_copy, (config, args)), Val(length(args)))...;
            kwargs...
        )

        primal = if Enzyme.EnzymeRules.needs_primal(config)
            res[1]
        else
            nothing
        end

        shadow = if Enzyme.EnzymeRules.needs_shadow(config)
            _make_solution_zero(res[1])::RT
        else
            nothing
        end
        tup = if Enzyme.EnzymeRules.needs_shadow(config)
            (shadow, res[2])
        else
            nothing
        end
        return Enzyme.EnzymeRules.augmented_rule_return_type(config, RTA)(primal, shadow, tup)
    end

    function Enzyme.EnzymeRules.reverse(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, ::Union{Type{Duplicated{RT}}, Type{MixedDuplicated{RT}}}, tape, prob,
            sensealg::Union{
                Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm},
            },
            u0, p, args...; kwargs...
        ) where {RT}

        if Enzyme.EnzymeRules.needs_shadow(config)
            dres, clos = tape
            dres = dres::RT
            dargs = clos(dres)
            # Mirror the `diff_tunables` choice the inner adjoint will make. When
            # the user passes a concrete sensealg, honor its `diff_tunables` field.
            # When the outer sensealg is `nothing` (default), `_concrete_solve_adjoint`
            # delegates to `automatic_sensealg_choice`, which picks
            # `diff_tunables = Val(false)` whenever `prob.p` is a SciMLStructure
            # with a non-empty `caches` field (e.g. an MTKParameters tied to an
            # SCCNonlinearProblem's `explicitfuns!` coupling). Reproducing that
            # predicate here lets the accumulator walk every non-Tunable field of
            # a structured `darg` so the meaningful cotangent isn't dropped.
            diff_tunables = let s = sensealg.val, pv = p.val
                if s isa DiffEqBase.AbstractSensitivityAlgorithm &&
                        hasproperty(s, :diff_tunables)
                    !(getproperty(s, :diff_tunables) isa Val{false})
                else
                    !(
                        SciMLStructures.isscimlstructure(pv) &&
                            !(pv isa AbstractArray) &&
                            hasfield(typeof(pv), :caches) &&
                            !isempty(pv.caches)
                    )
                end
            end
            for (darg, ptr) in zip(dargs, (func, prob, sensealg, u0, p, args...))
                if ptr isa Enzyme.Const
                    continue
                end
                if darg == ChainRulesCore.NoTangent()
                    continue
                end
                if ptr isa MixedDuplicated
                    _accum_tangent!(ptr.dval[], darg; diff_tunables)
                else
                    _accum_tangent!(ptr.dval, darg; diff_tunables)
                end
            end
            Enzyme.make_zero!(dres.u)
        end

        return ntuple(Returns(nothing), Val(length(args) + 4))
    end
end

end
