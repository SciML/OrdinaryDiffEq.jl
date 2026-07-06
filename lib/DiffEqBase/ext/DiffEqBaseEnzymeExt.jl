module DiffEqBaseEnzymeExt

@static if isempty(VERSION.prerelease)
    using DiffEqBase
    import DiffEqBase: value
    using Enzyme
    import Enzyme: Const, MixedDuplicated
    using ChainRulesCore


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

    # Note these following functions are generally not considered user facing from within Enzyme.
    # They enable additional performance/usability here (e.g. inactive kwargs).
    # Contact wsmoses@ before modifying (and beware their semantics may change without semver).

    Enzyme.EnzymeRules.inactive_kwarg(::typeof(DiffEqBase.solve_up), prob, sensalg::Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}, u0, p, args...; kwargs...) = nothing

    Enzyme.EnzymeRules.has_easy_rule(::typeof(DiffEqBase.solve_up), prob, sensalg::Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}, u0, p, args...; kwargs...) = nothing

    # `sensealg` is semantically inactive (a configuration token, not a value
    # whose tangent has meaning), but Enzyme's `set_runtime_activity(Reverse)`
    # can still promote module-level `const _MTK_SENSEALG = GaussAdjoint(...)`
    # references — or any sensealg captured in an `ODEProblem` field — to
    # `Duplicated{<:AbstractSensitivityAlgorithm}` before this rule is matched.
    # Accept any `Enzyme.Annotation` wrapping a sensealg value (Const /
    # Duplicated / MixedDuplicated / Active) and read only `sensealg.val`;
    # the shadow has no meaning. The inner-type constraint keeps
    # `SensitivityADPassThrough` (which is `<:AbstractDEAlgorithm`, not
    # `<:AbstractSensitivityAlgorithm`) off this rule so that case still
    # falls through to direct AD of `solve_up`.
    function Enzyme.EnzymeRules.augmented_primal(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, RTA::Union{Type{Duplicated{RT}}, Type{MixedDuplicated{RT}}}, prob,
            sensealg::Enzyme.Annotation{<:Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}},
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
            mz = Enzyme.make_zero(res[1])
            if Base.isabstracttype(RT) && (Enzyme.guess_activity(Core.Typeof(res[1]), Enzyme.Reverse) <: Enzyme.MixedDuplicated)
                Ref(mz)::RT
            else
                mz::RT
            end
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

    # Accumulate a reverse-pass cotangent `darg` into the Enzyme shadow `dval`.
    # Ordinary array shadows take a broadcasted add. An `MTKParameters` shadow
    # cannot be broadcast against its structural `NamedTuple` tangent, so its
    # differentiable buffers (`tunable`, `caches`) are accumulated directly.
    @inline _accumulate_tangent!(dval, darg) = (dval .+= darg; nothing)
    @inline _accumulate_tangent!(dval, ::Nothing) = nothing
    @inline function _accumulate_tangent!(dval, darg::NamedTuple)
        if hasproperty(darg, :tunable) && darg.tunable !== nothing
            dval.tunable .+= darg.tunable
        end
        if hasproperty(darg, :caches) && darg.caches !== nothing
            for (c, dc) in zip(dval.caches, darg.caches)
                dc === nothing || (c .+= dc)
            end
        end
        nothing
    end

    function Enzyme.EnzymeRules.reverse(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, ::Union{Type{Duplicated{RT}}, Type{MixedDuplicated{RT}}}, tape, prob,
            sensealg::Enzyme.Annotation{<:Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm}},
            u0, p, args...; kwargs...
        ) where {RT}

        if Enzyme.EnzymeRules.needs_shadow(config)
            dres, clos = tape
            dres = dres::RT
            # augmented_primal returns a Ref-wrapped shadow for the abstract
            # MixedDuplicated case; unwrap before handing the cotangent to the
            # ChainRules pullback (which expects the bare ODESolution tangent).
            dval = dres isa Base.RefValue ? dres[] : dres
            dargs = clos(dval)
            # `sensealg` is inactive (see augmented_primal note); skip its slot
            # whether it arrived as Const or as a runtime-activity-promoted
            # Duplicated/MixedDuplicated/Active.
            rta = Enzyme.EnzymeRules.runtime_activity(config) # detect runtime activity
            for (darg, ptr) in zip(dargs, (func, prob, sensealg, u0, p, args...))
                if ptr isa Enzyme.Const
                    continue
                end
                if ptr === sensealg
                    continue
                end
                if darg == ChainRulesCore.NoTangent()
                    continue
                end
                # Under `set_runtime_activity` (detected by rta), a runtime-inactive value arrives as
                # Duplicated/MixedDuplicated with its shadow ALIASING the primal
                # (`dval === val`). Accumulating the cotangent into such a shadow writes
                # gradient values into the caller's primal data (e.g. the `u0` of a
                # `Const` problem whose array was reused via `remake`), silently
                # corrupting subsequent calls. Skip them: a runtime-inactive value
                # accumulates nowhere, exactly as if it were `Const`.
                if ptr isa MixedDuplicated
                    rta && ptr.dval[] === ptr.val && continue
                    _accumulate_tangent!(ptr.dval[], darg)
                else
                    rta && ptr.dval === ptr.val && continue
                    _accumulate_tangent!(ptr.dval, darg)
                end
            end
            Enzyme.make_zero!(dval.u)
        end

        return ntuple(Returns(nothing), Val(length(args) + 4))
    end
end

end
